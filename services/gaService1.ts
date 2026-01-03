import { NetworkState, GAParameters, GAResult, BusNode, BranchLine, SelectionMethod, NodeType } from '../types';

// Complex number class
class Complex {
  constructor(public real: number, public imag: number) {}

  add(c: Complex): Complex {
    return new Complex(this.real + c.real, this.imag + c.imag);
  }

  sub(c: Complex): Complex {
    return new Complex(this.real - c.real, this.imag - c.imag);
  }

  mult(c: Complex): Complex {
    return new Complex(
      this.real * c.real - this.imag * c.imag,
      this.real * c.imag + this.imag * c.real
    );
  }

  magnitude(): number {
    return Math.sqrt(this.real * this.real + this.imag * this.imag);
  }

  angle(): number {
    return Math.atan2(this.imag, this.real);
  }
}

// Build Y-bus matrix from network
function buildYbus(network: NetworkState): Complex[][] {
  const n = network.nodes.length;
  const Y = Array(n).fill(0).map(() => Array(n).fill(0).map(() => new Complex(0, 0)));

  network.edges.forEach(edge => {
    const fromIdx = network.nodes.findIndex(node => node.id === edge.from);
    const toIdx = network.nodes.findIndex(node => node.id === edge.to);

    if (fromIdx === -1 || toIdx === -1) return;

    const y = new Complex(
      edge.resistance / (edge.resistance * edge.resistance + edge.reactance * edge.reactance),
      -edge.reactance / (edge.resistance * edge.resistance + edge.reactance * edge.reactance)
    );
    const b = new Complex(0, edge.susceptance);

    Y[fromIdx][fromIdx] = Y[fromIdx][fromIdx].add(y).add(b);
    Y[toIdx][toIdx] = Y[toIdx][toIdx].add(y).add(b);
    Y[fromIdx][toIdx] = Y[fromIdx][toIdx].sub(y);
    Y[toIdx][fromIdx] = Y[toIdx][fromIdx].sub(y);
  });

  return Y;
}

// Calculate power at a bus
function calculatePower(V: number[], theta: number[], Ybus: Complex[][], busIdx: number): { P: number; Q: number } {
  let P = 0, Q = 0;
  const n = V.length;

  for (let j = 0; j < n; j++) {
    const Vj = V[j];
    const thetaIJ = theta[busIdx] - theta[j];
    const Y = Ybus[busIdx][j];
    const G = Y.real;
    const B = Y.imag;

    P += Vj * (G * Math.cos(thetaIJ) + B * Math.sin(thetaIJ));
    Q += Vj * (G * Math.sin(thetaIJ) - B * Math.cos(thetaIJ));
  }

  return { P: V[busIdx] * P, Q: V[busIdx] * Q };
}

// Fitness function for GA
function fitness(chromosome: number[], network: NetworkState, Ybus: Complex[][], vMin: number, vMax: number, alpha: number, beta: number, gamma: number): number {
  const n = network.nodes.length;
  const V = chromosome.slice(0, n);
  const theta = chromosome.slice(n);

  let penalty = 0;

  network.nodes.forEach((bus, i) => {
    if (bus.type === NodeType.SLACK) return;

    const calc = calculatePower(V, theta, Ybus, i);
    const Pspec = bus.pGen - bus.pLoad;
    const Qspec = bus.qGen - bus.qLoad;

    // Power balance penalties
    penalty += alpha * Math.pow(calc.P - Pspec, 2); // delta P violations
    if (bus.type === NodeType.PQ) {
      penalty += beta * Math.pow(calc.Q - Qspec, 2); // delta Q violations
    }

    // Voltage constraint penalty
    if (V[i] < vMin || V[i] > vMax) {
      penalty += 1000 * Math.pow(Math.max(vMin - V[i], V[i] - vMax, 0), 2);
    }

    // Reactive power generation (Qg) constraint violations for generators
    if (bus.type === NodeType.PV || bus.type === NodeType.SLACK) {
      // For generators, Qg should be within reasonable limits
      // This is a simplified constraint - in practice you'd have specific limits
      const qgViolation = Math.max(0, Math.abs(calc.Q) - 2.0); // Assume |Qg| <= 2.0 p.u. as example
      penalty += gamma * Math.pow(qgViolation, 2);
    }
  });

  return penalty;
}

// Tournament selection
function tournament(population: number[][], fitnessScores: number[], k: number): number[] {
  let best = null;
  let bestFitness = Infinity;

  for (let i = 0; i < k; i++) {
    const idx = Math.floor(Math.random() * population.length);
    if (fitnessScores[idx] < bestFitness) {
      bestFitness = fitnessScores[idx];
      best = population[idx];
    }
  }

  return best!;
}

// Roulette wheel selection
function rouletteWheel(population: number[][], fitnessScores: number[]): number[] {
  // Convert fitness to selection probabilities (lower fitness = higher probability)
  const minFitness = Math.min(...fitnessScores);
  const maxFitness = Math.max(...fitnessScores);
  const range = maxFitness - minFitness || 1; // Avoid division by zero

  // Invert fitness values (higher fitness becomes lower probability)
  const invertedFitness = fitnessScores.map(f => maxFitness - f + minFitness);
  const totalFitness = invertedFitness.reduce((sum, f) => sum + f, 0);

  const probabilities = invertedFitness.map(f => f / totalFitness);
  const cumulative = probabilities.reduce((acc, p, i) => {
    acc.push((acc[i - 1] || 0) + p);
    return acc;
  }, [] as number[]);

  const rand = Math.random();
  const selectedIndex = cumulative.findIndex(cum => cum >= rand);
  return population[selectedIndex !== -1 ? selectedIndex : 0];
}

// Rank-based selection
function rankBased(population: number[][], fitnessScores: number[]): number[] {
  // Sort by fitness (ascending - lower fitness is better)
  const sortedIndices = fitnessScores
    .map((f, i) => ({ fitness: f, index: i }))
    .sort((a, b) => a.fitness - b.fitness)
    .map(item => item.index);

  // Assign ranks (lower rank number = better fitness)
  const ranks = new Array(population.length);
  sortedIndices.forEach((originalIndex, rank) => {
    ranks[originalIndex] = rank + 1; // 1-based rank
  });

  // Calculate selection probabilities based on ranks
  const totalRank = (population.length * (population.length + 1)) / 2;
  const probabilities = ranks.map(rank => (population.length - rank + 1) / totalRank);

  const cumulative = probabilities.reduce((acc, p, i) => {
    acc.push((acc[i - 1] || 0) + p);
    return acc;
  }, [] as number[]);

  const rand = Math.random();
  const selectedIndex = cumulative.findIndex(cum => cum >= rand);
  return population[selectedIndex !== -1 ? selectedIndex : 0];
}

// Stochastic universal sampling
function stochasticUniversal(population: number[][], fitnessScores: number[]): number[][] {
  // This returns multiple parents for SUS
  const numParents = 2; // Typically select 2 parents

  // Convert fitness to selection probabilities (same as roulette wheel)
  const minFitness = Math.min(...fitnessScores);
  const maxFitness = Math.max(...fitnessScores);
  const range = maxFitness - minFitness || 1;

  const invertedFitness = fitnessScores.map(f => maxFitness - f + minFitness);
  const totalFitness = invertedFitness.reduce((sum, f) => sum + f, 0);

  const probabilities = invertedFitness.map(f => f / totalFitness);
  const cumulative = probabilities.reduce((acc, p, i) => {
    acc.push((acc[i - 1] || 0) + p);
    return acc;
  }, [] as number[]);

  const selectedParents: number[][] = [];
  const spacing = 1.0 / numParents;

  for (let i = 0; i < numParents; i++) {
    const rand = (Math.random() + i) * spacing;
    const selectedIndex = cumulative.findIndex(cum => cum >= rand);
    const actualIndex = selectedIndex !== -1 ? selectedIndex : 0;
    selectedParents.push(population[actualIndex]);
  }

  return selectedParents;
}

// Unified selection function
function selectParents(population: number[][], fitnessScores: number[], method: SelectionMethod, tournamentSize: number = 5): number[][] {
  switch (method) {
    case SelectionMethod.TOURNAMENT:
      return [tournament(population, fitnessScores, tournamentSize), tournament(population, fitnessScores, tournamentSize)];

    case SelectionMethod.ROULETTE_WHEEL:
      return [rouletteWheel(population, fitnessScores), rouletteWheel(population, fitnessScores)];

    case SelectionMethod.RANK_BASED:
      return [rankBased(population, fitnessScores), rankBased(population, fitnessScores)];

    case SelectionMethod.STOCHASTIC_UNIVERSAL:
      return stochasticUniversal(population, fitnessScores);

    default:
      return [tournament(population, fitnessScores, 5), tournament(population, fitnessScores, 5)];
  }
}

// Mutation
function mutate(chromosome: number[], rate: number, vMin: number, vMax: number, network: NetworkState): number[] {
  const n = network.nodes.length;
  const mutated = [...chromosome];

  for (let i = 0; i < chromosome.length; i++) {
    if (Math.random() < rate) {
      if (i < n) {
        // Voltage mutation
        if (network.nodes[i].type !== NodeType.SLACK) {
          mutated[i] += (Math.random() - 0.5) * 0.05;
          mutated[i] = Math.max(vMin, Math.min(vMax, mutated[i]));
        }
      } else {
        // Angle mutation
        const busIdx = i - n;
        if (network.nodes[busIdx].type !== NodeType.SLACK) {
          mutated[i] += (Math.random() - 0.5) * 0.1;
        }
      }
    }
  }

  return mutated;
}

// Main GA function
export async function calculateLoadFlowGA(
  network: NetworkState,
  params: GAParameters,
  progressCallback?: (generation: number, fitness: number) => void
): Promise<GAResult> {
  const { populationSize, maxGenerations, crossoverRate, mutationRate, convergenceThreshold, minVoltage, maxVoltage, selectionMethod, tournamentSize } = params;
  const n = network.nodes.length;
  const Ybus = buildYbus(network);

  // Initialize population
  let population: number[][] = [];
  for (let i = 0; i < populationSize; i++) {
    const chromosome: number[] = [];
    // Voltage magnitudes
    for (let j = 0; j < n; j++) {
      chromosome.push(network.nodes[j].type === NodeType.SLACK ?
                     network.nodes[j].voltage :
                     (minVoltage + Math.random() * (maxVoltage - minVoltage)));
    }
    // Angles (radians)
    for (let j = 0; j < n; j++) {
      chromosome.push(network.nodes[j].type === NodeType.SLACK ? 0 :
                     (Math.random() * 0.2 - 0.1));
    }
    population.push(chromosome);
  }

  const fitnessHistory: number[] = [];
  let bestSolution = null;
  let bestFitness = Infinity;

  for (let gen = 0; gen < maxGenerations; gen++) {
    // Evaluate fitness
    const fitnessScores = population.map(ind => fitness(ind, network, Ybus, minVoltage, maxVoltage, params.alpha, params.beta, params.gamma));

    const minFitness = Math.min(...fitnessScores);
    const minIdx = fitnessScores.indexOf(minFitness);

    if (minFitness < bestFitness) {
      bestFitness = minFitness;
      bestSolution = [...population[minIdx]];
    }

    fitnessHistory.push(bestFitness);

    if (progressCallback) {
      progressCallback(gen + 1, bestFitness);
    }

    // Check convergence
    if (bestFitness < convergenceThreshold) break;

    // Selection and create new population
    const newPop: number[][] = [];
    const eliteCount = Math.floor(populationSize * 0.1);

    // Elitism
    const sortedIndices = fitnessScores
      .map((f, i) => ({ f, i }))
      .sort((a, b) => a.f - b.f)
      .slice(0, eliteCount)
      .map(x => x.i);

    sortedIndices.forEach(idx => newPop.push([...population[idx]]));

    // Generate rest of population using selected method
    while (newPop.length < populationSize) {
      const parents = selectParents(population, fitnessScores, selectionMethod, tournamentSize || 5);
      const parent1 = parents[0];
      const parent2 = parents[1];

      let child1 = [...parent1];
      let child2 = [...parent2];

      // Crossover
      if (Math.random() < crossoverRate) {
        const point = Math.floor(Math.random() * parent1.length);
        child1 = [...parent1.slice(0, point), ...parent2.slice(point)];
        child2 = [...parent2.slice(0, point), ...parent1.slice(point)];
      }

      // Mutation
      child1 = mutate(child1, mutationRate, minVoltage, maxVoltage, network);
      child2 = mutate(child2, mutationRate, minVoltage, maxVoltage, network);

      newPop.push(child1);
      if (newPop.length < populationSize) newPop.push(child2);
    }

    population = newPop;

    // Allow UI to update
    if (gen % 10 === 0) {
      await new Promise(resolve => setTimeout(resolve, 0));
    }
  }

  return {
    solution: {
      voltages: bestSolution!.slice(0, n),
      angles: bestSolution!.slice(n).map(a => a * 180 / Math.PI) // Convert to degrees
    },
    fitness: bestFitness,
    generations: fitnessHistory.length,
    fitnessHistory,
    converged: bestFitness < convergenceThreshold
  };
}
