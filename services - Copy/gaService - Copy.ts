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

function buildYbus(network: NetworkState): Complex[][] {
  const n = network.nodes.length;
  const Y = Array(n).fill(0).map(() => Array(n).fill(0).map(() => new Complex(0, 0)));

  network.edges.forEach(edge => {
    const fromIdx = network.nodes.findIndex(node => node.id === edge.from);
    const toIdx = network.nodes.findIndex(node => node.id === edge.to);

    if (fromIdx === -1 || toIdx === -1) return;

    const Z_mag_sq = edge.resistance * edge.resistance + edge.reactance * edge.reactance;
    
    // FIX: Skip edge entirely if impedance is zero
    if (Z_mag_sq === 0 || Z_mag_sq < 1e-10) {
      console.warn(`Zero or near-zero impedance for edge ${edge.from}-${edge.to}, skipping`);
      return;  // Early exit - don't process this edge
    }
    
    const y = new Complex(
      edge.resistance / Z_mag_sq,
      -edge.reactance / Z_mag_sq
    );
    
    const y_shunt = new Complex(0, edge.susceptance / 2);

    Y[fromIdx][fromIdx] = Y[fromIdx][fromIdx].add(y).add(y_shunt);
    Y[toIdx][toIdx] = Y[toIdx][toIdx].add(y).add(y_shunt);
    Y[fromIdx][toIdx] = Y[fromIdx][toIdx].sub(y);
    Y[toIdx][fromIdx] = Y[toIdx][fromIdx].sub(y);
  });

  return Y;
}

// Calculate P and Q injections at all buses
function calculatePQ(V: number[], theta: number[], Ybus: Complex[][]): { P: number[], Q: number[] } {
  const n = V.length;
  const P = new Array(n).fill(0);
  const Q = new Array(n).fill(0);

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const G = Ybus[i][j].real;
      const B = Ybus[i][j].imag;
      const thetaIJ = theta[i] - theta[j];
      
      P[i] += V[i] * V[j] * (G * Math.cos(thetaIJ) + B * Math.sin(thetaIJ));
      Q[i] += V[i] * V[j] * (G * Math.sin(thetaIJ) - B * Math.cos(thetaIJ));
    }
  }

  return { P , Q };
}


// Calculate total system losses
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

// Fitness function with correct formulation
function fitness(
  chromosome: number[], 
  network: NetworkState, 
  Ybus: Complex[][], 
  vMin: number, 
  vMax: number,
  alpha: number,
  beta: number
): number {
  const n = network.nodes.length;
  let idx = 0;

  // Decode chromosome
  const V = new Array(n);
  const theta = new Array(n);

  // Slack bus (first bus by convention) - fixed values
  V[0] = network.nodes[0].voltage;
  theta[0] = 0;

  // PV buses - voltage fixed, angle optimized
  network.nodes.forEach((node, i) => {
    if (i === 0) return; // Skip slack
    
    if (node.type === NodeType.PV) {
      V[i] = node.voltage; // Fixed voltage
      theta[i] = chromosome[idx++]; // Optimize angle
    }
  });

  // PQ buses - both voltage and angle optimized
  network.nodes.forEach((node, i) => {
    if (i === 0) return; // Skip slack
    
    if (node.type === NodeType.PQ) {
      V[i] = chromosome[idx++]; // Optimize voltage
      theta[i] = chromosome[idx++]; // Optimize angle
    }
  });

  // Calculate actual P and Q at each bus

  const { P: Pcalc, Q: Qcalc } = calculatePQ(V, theta, Ybus);

  // Calculate penalty for power balance violations
  let penalty = 0;

  network.nodes.forEach((node, i) => {
    // Net injection = Generation - Load
    const Pnet_spec = (node.pGen - node.pLoad)/100;
    const Qnet_spec = (node.qGen - node.qLoad)/100;

    if (node.type === NodeType.SLACK) {
      // Slack bus: Q mismatch only (P is free to balance system)
      //const Qmismatch = Qcalc[i] - Qnet_spec;
      //penalty += beta * Qmismatch * Qmismatch;
      
    } else if (node.type === NodeType.PV) {
      // PV bus: P must match, Q is free (within limits)
      const Pmismatch = Pcalc[i] - Pnet_spec;
      //console.log(`Bus ${node.pLoad} ${Pcalc[i]} ${Pcalc[i]+node.pLoad }P mismatch: ${Pmismatch} , Q mismatch: ${Pnet_spec}`);
      penalty += Math.pow(Pmismatch, 2)  ;
      
      // Soft Q limits (typical generator limits)
      //const Qgen = Qcalc[i] + node.qLoad; // Generator reactive power
      //if (Qgen < -2.0) {
       // penalty += 1000 * Math.pow(Qgen + 2.0, 2);
      //} else if (Qgen > 2.0) {
        //penalty += 1000 * Math.pow(Qgen - 2.0, 2);
      //}
      
    } else { // PQ bus
      // PQ bus: both P and Q must match
      const Pmismatch = Pcalc[i] - Pnet_spec;
      const Qmismatch = Qcalc[i] - Qnet_spec;
      penalty += Math.pow(Pmismatch, 2);
      penalty += Math.pow(Qmismatch, 2);
    }

    // Voltage limits for PQ buses
    if (node.type === NodeType.PQ) {
      // Voltage constraint penalty
                if (V[i] < vMin || V[i] > vMax) {
                    penalty += 1000 * Math.pow(Math.max(vMin - V[i], V[i] - vMax, 0), 2);
                }
    }
  });

  // Add losses as part of objective (minimize losses + constraint violations)
  const losses = calculateLosses(V, theta, network);
  console.log(`Losses: ${losses}`);
  return penalty + losses;
}

// Gaussian random number generator
function gaussianRandom(): number {
  // FIX: Prevent u1 or u2 from being exactly 0
  const u1 = Math.random()   ;        
  const u2 = Math.random() ;
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}

// Tournament selection
function tournamentSelection(population: number[][], fitnessScores: number[], k?: number): number[] {
  // FIX: Adaptive tournament size based on population
  if (!k) {
    k = Math.max(2, Math.min(5, Math.floor(population.length / 20)));
  }
  
  let bestIdx = Math.floor(Math.random() * population.length);
  let bestFitness = fitnessScores[bestIdx];

  for (let i = 1; i < k; i++) {
    const idx = Math.floor(Math.random() * population.length);
    if (fitnessScores[idx] < bestFitness) {
      bestFitness = fitnessScores[idx];
      bestIdx = idx;
    }
  }

  return [...population[bestIdx]];
}

// Crossover (arithmetic/blend)
function crossover(parent1: number[], parent2: number[], rate: number): [number[], number[]] {
  if (Math.random() > rate) {
    return [[...parent1], [...parent2]];
  }

  const alpha = 0.5;
  const child1 = parent1.map((val, i) => alpha * val + (1 - alpha) * parent2[i]);
  const child2 = parent1.map((val, i) => (1 - alpha) * val + alpha * parent2[i]);

  return [child1, child2];
}

// Mutation
function mutate(chromosome: number[], rate: number, network: NetworkState, vMin: number, vMax: number): number[] {
  const mutated = [...chromosome];
  let idx = 0;

  const angleSigma = 0.3;
  const voltageSigma = 0.02;

  // PV buses - only angles
  network.nodes.forEach((node, i) => {
    if (i === 0) return;
    
    if (node.type === NodeType.PV) {
      if (Math.random() < rate) {
        mutated[idx] +=  (Math.random() - 0.5) * 0.1;
      }
      idx++;
    }
  });

  // PQ buses - voltages and angles
  network.nodes.forEach((node, i) => {
    if (i === 0) return;
    
    if (node.type === NodeType.PQ) {
      // Voltage mutation
      if (Math.random() < rate) {
        mutated[idx] += (Math.random() - 0.5) * 0.05;
        mutated[idx] = Math.max(vMin, Math.min(vMax, mutated[idx]));
      }
      idx++;
      
      // Angle mutation
      if (Math.random() < rate) {
        mutated[idx] +=(Math.random() - 0.5) * 0.1;
      }
      idx++;
    }
  });

  return mutated;
}

// Main GA function
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
    alpha,
    beta
  } = params;

  // Build Y-bus
  const Ybus = buildYbus(network);
  // Determine chromosome size
  let chromSize = 0;
  network.nodes.forEach((node, i) => {
    if (i === 0) return; // Skip slack
    if (node.type === NodeType.PV) {
      chromSize += 1; // angle only
    } else if (node.type === NodeType.PQ) {
      chromSize += 2; // voltage and angle
    }
  });

  // Initialize population
  let population: number[][] = [];
// FIX: Build chromosome in single pass with clear indexing
for (let i = 0; i < populationSize; i++) {
  const chromosome: number[] = [];

  network.nodes.forEach((node, j) => {
    if (j === 0) return;  // Skip slack bus
    
    if (node.type === NodeType.PV) {
      // PV bus: only angle
      chromosome.push((Math.random() - 0.5) * 1.0);
    } else if (node.type === NodeType.PQ) {
      // PQ bus: voltage then angle
      chromosome.push(minVoltage + Math.random() * (maxVoltage - minVoltage));
      chromosome.push((Math.random() - 0.5) * 1.0);
    }
  });

  population.push(chromosome);
}

  const fitnessHistory: number[] = [];
  let bestSolution: number[] | null = null;
  let bestFitness = Infinity;
  let bestVoltages: number[] = [];
  let bestAngles: number[] = [];

  for (let gen = 0; gen < maxGenerations; gen++) {
    // Evaluate fitness
    const fitnessScores = population.map(ind => 
      fitness(ind, network, Ybus, minVoltage, maxVoltage, alpha, beta)
    );

    const minFitness = Math.min(...fitnessScores);
    const minIdx = fitnessScores.indexOf(minFitness);

    if (minFitness < bestFitness) {
      bestFitness = minFitness;
      bestSolution = [...population[minIdx]];

      // Decode best solution for output
      const n = network.nodes.length;
      const V = new Array(n);
      const theta = new Array(n);
      let idx = 0;

      V[0] = network.nodes[0].voltage;
      theta[0] = 0;

      network.nodes.forEach((node, i) => {
        if (i === 0) return;
        if (node.type === NodeType.PV) {
          V[i] = node.voltage;
          theta[i] = bestSolution![idx++];
        }
      });

      network.nodes.forEach((node, i) => {
        if (i === 0) return;
        if (node.type === NodeType.PQ) {
          V[i] = bestSolution![idx++];
          theta[i] = bestSolution![idx++];
        }
      });

      bestVoltages = V;
      bestAngles = theta;
    }

    fitnessHistory.push(bestFitness);

    if (progressCallback) {
      progressCallback(gen + 1, bestFitness);
    }

    // Check convergence
    if (bestFitness < convergenceThreshold) break;

    // Create new population
    const newPopulation: number[][] = [];

    // Elitism - keep best solution
    newPopulation.push([...population[minIdx]]);

    // Generate rest of population
    while (newPopulation.length < populationSize) {
      const parent1 = tournamentSelection(population, fitnessScores);
      const parent2 = tournamentSelection(population, fitnessScores);

      let [child1, child2] = crossover(parent1, parent2, crossoverRate);

      child1 = mutate(child1, mutationRate, network, minVoltage, maxVoltage);
      child2 = mutate(child2, mutationRate, network, minVoltage, maxVoltage);

      newPopulation.push(child1);
      if (newPopulation.length < populationSize) {
        newPopulation.push(child2);
      }
    }

    population = newPopulation;

    // Allow UI updates
    if (gen % 10 === 0) {
      await new Promise(resolve => setTimeout(resolve, 0));
    }
  }

  // Calculate final P and Q values
  const { P: Pcalc, Q: Qcalc } = calculatePQ(bestVoltages, bestAngles, Ybus);

  return {
    solution: {
      voltages: bestVoltages,
      angles: bestAngles.map(a => a * 180 / Math.PI),
      generations: {
        pGen: Pcalc.map((p, i) => p + network.nodes[i].pLoad),
        qGen: Qcalc.map((q, i) => q + network.nodes[i].qLoad)
      }
    },
    fitness: bestFitness,
    generations: fitnessHistory.length,
    fitnessHistory,
    converged: bestFitness < convergenceThreshold
  };
}
