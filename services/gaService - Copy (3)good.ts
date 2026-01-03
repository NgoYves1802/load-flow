import { NetworkState, GAParameters, GAResult, NodeType, SelectionMethod } from '../types';

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

// Build Y-bus matrix - EXACT HTML logic
function buildYbus(network: NetworkState): Complex[][] {
  const n = network.nodes.length;
  const Y = Array(n).fill(0).map(() => Array(n).fill(0).map(() => new Complex(0, 0)));

  network.edges.forEach(edge => {
    const fromIdx = network.nodes.findIndex(node => node.id === edge.from);
    const toIdx = network.nodes.findIndex(node => node.id === edge.to);

    if (fromIdx === -1 || toIdx === -1) return;

    // EXACT HTML calculation
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

// Calculate power at a bus - EXACT HTML logic
function calculatePower(V: number[], theta: number[], busIdx: number, Ybus: Complex[][]): { P: number, Q: number } {
  let P = 0;
  let Q = 0;
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
function calculateLosses(V: number[], theta: number[], network: NetworkState): number {
  let totalLoss = 0;
  
  network.edges.forEach(edge => {
    const fromIdx = network.nodes.findIndex(node => node.id === edge.from);
    const toIdx = network.nodes.findIndex(node => node.id === edge.to);
    
    if (fromIdx === -1 || toIdx === -1) return;
    
    const Z_mag_sq = edge.resistance * edge.resistance + edge.reactance * edge.reactance;
    if (Z_mag_sq === 0) return;
    
    const G = edge.resistance / Z_mag_sq;
    const thetaIJ = theta[fromIdx] - theta[toIdx];
    
    totalLoss += G * (
      V[fromIdx] * V[fromIdx] + 
      V[toIdx] * V[toIdx] - 
      2 * V[fromIdx] * V[toIdx] * Math.cos(thetaIJ)
    );
  });
  
  return totalLoss;
}
// Fitness function - with baseMVA conversion
function fitness(chromosome: number[], network: NetworkState, Ybus: Complex[][], vMin: number, vMax: number, baseMVA: number, alpha: number, beta: number, gamma: number): number {
  const n = network.nodes.length;
  const V = [...chromosome.slice(0, n)];
  const theta = [...chromosome.slice(n)];

  let error = 0;
  let penalty = 0;

  network.nodes.forEach((bus, i) => {
    if (bus.type === NodeType.SLACK) return;

    const calc = calculatePower(V, theta, i, Ybus);
    
    // Convert MW/MVAr to p.u. for comparison
    const Pspec = (bus.pGen - bus.pLoad) / baseMVA;
    const Qspec = (bus.qGen - bus.qLoad) / baseMVA;

    error += alpha * Math.pow(calc.P - Pspec, 2);
    if (bus.type === NodeType.PQ) {
      error += beta  * Math.pow(calc.Q - Qspec, 2);
    }

    // Voltage constraint penalty
    if (V[i] < vMin || V[i] > vMax) {
      penalty += 1000 * Math.pow(Math.max(vMin - V[i], V[i] - vMax, 0), 2);
    }
  });
 const losses = calculateLosses(V, theta, network);
  console.log(`Losses: ${losses}`);
  return error + penalty + gamma * losses;
}

// Tournament selection - EXACT HTML logic
function tournament(population: number[][], fitnessScores: number[], k: number): number[] {
  let best: number[] | null = null;
  let bestFitness = Infinity;

  for (let i = 0; i < k; i++) {
    const idx = Math.floor(Math.random() * population.length);
    if (fitnessScores[idx] < bestFitness) {
      bestFitness = fitnessScores[idx];
      best = population[idx];
    }
  }

  return [...best!];
}

// Roulette wheel selection
function rouletteWheel(population: number[][], fitnessScores: number[]): number[] {
  // Fitness scores are minimization problems, so invert them
  const maxFitness = Math.max(...fitnessScores);
  const invertedFitness = fitnessScores.map(f => maxFitness - f + 1e-6); // Add small epsilon to avoid zero
  const totalFitness = invertedFitness.reduce((sum, f) => sum + f, 0);

  let pick = Math.random() * totalFitness;
  let currentSum = 0;

  for (let i = 0; i < population.length; i++) {
    currentSum += invertedFitness[i];
    if (currentSum >= pick) {
      return [...population[i]];
    }
  }

  // Fallback (should not reach here normally)
  return [...population[population.length - 1]];
}

// Rank-based selection
function rankBased(population: number[][], fitnessScores: number[]): number[] {
  // Sort indices by fitness (ascending for minimization)
  const sortedIndices = fitnessScores
    .map((f, i) => ({ f, i }))
    .sort((a, b) => a.f - b.f)
    .map(item => item.i);

  // Assign ranks (1 to N, where 1 is best)
  const ranks = new Array(population.length);
  sortedIndices.forEach((idx, rank) => {
    ranks[idx] = rank + 1;
  });

  // Calculate selection probabilities (linear ranking)
  const totalRankSum = ranks.reduce((sum, rank) => sum + rank, 0);
  const probabilities = ranks.map(rank => rank / totalRankSum);

  // Cumulative probabilities
  const cumulative = probabilities.reduce((acc, prob, i) => {
    acc.push((acc[i - 1] || 0) + prob);
    return acc;
  }, [] as number[]);

  let pick = Math.random();
  for (let i = 0; i < cumulative.length; i++) {
    if (pick <= cumulative[i]) {
      return [...population[i]];
    }
  }

  // Fallback
  return [...population[population.length - 1]];
}

// Stochastic universal sampling
function stochasticUniversal(population: number[][], fitnessScores: number[]): number[] {
  // Fitness scores are minimization problems, so invert them
  const maxFitness = Math.max(...fitnessScores);
  const invertedFitness = fitnessScores.map(f => maxFitness - f + 1e-6);
  const totalFitness = invertedFitness.reduce((sum, f) => sum + f, 0);

  // Start at random point
  const start = Math.random() * (totalFitness / population.length);
  const pointers = [];
  for (let i = 0; i < population.length; i++) {
    pointers.push(start + i * (totalFitness / population.length));
  }

  let currentSum = 0;
  let pointerIndex = 0;

  for (let i = 0; i < population.length; i++) {
    currentSum += invertedFitness[i];
    while (pointerIndex < pointers.length && pointers[pointerIndex] <= currentSum) {
      return [...population[i]];
    }
  }

  // Fallback
  return [...population[population.length - 1]];
}

// Generic selection function that dispatches to the appropriate method
function selectParent(population: number[][], fitnessScores: number[], selectionMethod: SelectionMethod, tournamentSize?: number): number[] {
  switch (selectionMethod) {
    case SelectionMethod.TOURNAMENT:
      return tournament(population, fitnessScores, tournamentSize || 5);
    case SelectionMethod.ROULETTE_WHEEL:
      return rouletteWheel(population, fitnessScores);
    case SelectionMethod.RANK_BASED:
      return rankBased(population, fitnessScores);
    case SelectionMethod.STOCHASTIC_UNIVERSAL:
      return stochasticUniversal(population, fitnessScores);
    default:
      // Default to tournament if unknown method
      return tournament(population, fitnessScores, tournamentSize || 5);
  }
}
function gaussianRandom(): number {
  const u1 = Math.random() || 1e-10;
  const u2 = Math.random() || 1e-10;
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}
// Mutation - EXACT HTML logic
function mutate(chromosome: number[], rate: number, network: NetworkState, vMin: number, vMax: number): number[] {
  const n = network.nodes.length;
  const mutated = [...chromosome];

  for (let i = 0; i < chromosome.length; i++) {
    if (Math.random() < rate) {
      if (i < n) {
        // Voltage mutation
        if (network.nodes[i].type !== NodeType.SLACK) {
          mutated[i] += (Math.random() - 0.5) * 0.1;
          mutated[i] = Math.max(vMin, Math.min(vMax, mutated[i]));
        }
      } else {
        // Angle mutation
        const busIdx = i - n;
        if (network.nodes[busIdx].type !== NodeType.SLACK) {
          mutated[i] += gaussianRandom() * 0.2;
        }
      }
    }
  }

  return mutated;
}

// Main GA function - with baseMVA
export async function calculateLoadFlowGA(
  network: NetworkState,
  params: GAParameters,
  progressCallback?: (generation: number, fitness: number) => void
): Promise<GAResult> {
  const {
    populationSize,
    maxGenerations,
    crossoverRate,
    mutationRate,
    convergenceThreshold,
    minVoltage,
    maxVoltage,
    selectionMethod,
    tournamentSize,
    alpha,
    beta,
    gamma,
  } = params;

  const n = network.nodes.length;
  const baseMVA =  100; // Base power in MVA

  // Build Y-bus
  const Ybus = buildYbus(network);

  // Initialize population - EXACT HTML logic
  let population: number[][] = [];
  for (let i = 0; i < populationSize; i++) {
    const chromosome: number[] = [];
    
    // Voltage magnitudes
    for (let j = 0; j < n; j++) {
      chromosome.push(
        network.nodes[j].type === NodeType.SLACK 
          ? network.nodes[j].voltage 
          : 1.0
      );
    }
    
    // Angles (radians)
    for (let j = 0; j < n; j++) {
      chromosome.push(
        network.nodes[j].type === NodeType.SLACK 
          ? 0 
          : ((Math.random() - 0.5) * 1.0)
      );
    }
    
    population.push(chromosome);
  }

  const fitnessHistory: number[] = [];
  let bestSolution: number[] | null = null;
  let bestFitness = Infinity;

  for (let gen = 0; gen < maxGenerations; gen++) {
    // Evaluate fitness
    const fitnessScores = population.map(ind => fitness(ind, network, Ybus, minVoltage, maxVoltage, baseMVA, alpha, beta, gamma));

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

    if (bestFitness < convergenceThreshold) break;

    // Selection - EXACT HTML logic
    const newPop: number[][] = [];
    const eliteCount = Math.floor(populationSize * 0.1);

    // Elitism
    const sortedIndices = fitnessScores
      .map((f, i) => ({ f, i }))
      .sort((a, b) => a.f - b.f)
      .slice(0, eliteCount)
      .map(x => x.i);

    sortedIndices.forEach(idx => newPop.push([...population[idx]]));

    // Generate rest of population
    while (newPop.length < populationSize) {
      const parent1 = selectParent(population, fitnessScores, selectionMethod, tournamentSize);
      const parent2 = selectParent(population, fitnessScores, selectionMethod, tournamentSize);

      let child1 = [...parent1];
      let child2 = [...parent2];

      // Crossover
      if (Math.random() < crossoverRate) {
        const point = Math.floor(Math.random() * parent1.length);
        child1 = [...parent1.slice(0, point), ...parent2.slice(point)];
        child2 = [...parent2.slice(0, point), ...parent1.slice(point)];
      }

      // Mutation
      // Decrease mutation as generations progress
const adaptiveMutRate = mutationRate * (1 + Math.exp(-gen / 50));
child1 = mutate(child1, adaptiveMutRate, network, minVoltage, maxVoltage);
child2 = mutate(child2, adaptiveMutRate, network, minVoltage, maxVoltage);

      newPop.push(child1);
      if (newPop.length < populationSize) newPop.push(child2);
    }

    population = newPop;

    // Allow UI updates
    if (gen % 10 === 0) {
      await new Promise(resolve => setTimeout(resolve, 0));
    }
  }

  // Extract solution
  const finalV = bestSolution!.slice(0, n);
  const finalTheta = bestSolution!.slice(n);

  // Calculate final powers (in p.u.)
  const finalPowers = network.nodes.map((_, i) => 
    calculatePower(finalV, finalTheta, i, Ybus)
  );

  return {
    solution: {
      voltages: finalV,
      angles: finalTheta.map(a => a * 180 / Math.PI),
      generations: {
        // Convert p.u. back to MW/MVAr
        pGen: finalPowers.map((p, i) => (p.P * baseMVA) + network.nodes[i].pLoad),
        qGen: finalPowers.map((p, i) => (p.Q * baseMVA) + network.nodes[i].qLoad)
      }
    },
    fitness: bestFitness,
    generations: fitnessHistory.length,
    fitnessHistory,
    converged: bestFitness < convergenceThreshold
  };
}
