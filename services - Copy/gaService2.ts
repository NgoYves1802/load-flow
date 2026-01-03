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

    // FIXED: Correct calculation of y = 1/(R + jX)
    const Z_mag_sq = edge.resistance * edge.resistance + edge.reactance * edge.reactance;
    
    if (Z_mag_sq === 0) {
      console.warn(`Zero impedance for edge ${edge.from}-${edge.to}`);
      return;
    }
    
    const y = new Complex(
      edge.resistance / Z_mag_sq,   // G = R/(R²+X²)
      -edge.reactance / Z_mag_sq    // B = -X/(R²+X²)
    );
    
    // FIXED: Shunt susceptance split between both ends
    const y_shunt = new Complex(0, edge.susceptance / 2);

    Y[fromIdx][fromIdx] = Y[fromIdx][fromIdx].add(y).add(y_shunt);
    Y[toIdx][toIdx] = Y[toIdx][toIdx].add(y).add(y_shunt);
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

// ADDED: Calculate total system losses
function calculateTotalLosses(V: number[], theta: number[], network: NetworkState): number {
  let totalLoss = 0;
  
  network.edges.forEach(edge => {
    const fromIdx = network.nodes.findIndex(node => node.id === edge.from);
    const toIdx = network.nodes.findIndex(node => node.id === edge.to);
    
    if (fromIdx === -1 || toIdx === -1) return;
    
    const Z_mag_sq = edge.resistance * edge.resistance + edge.reactance * edge.reactance;
    if (Z_mag_sq === 0) return;
    
    const G = edge.resistance / Z_mag_sq;
    const thetaIJ = theta[fromIdx] - theta[toIdx];
    
    // Line loss = G * (Vi² + Vj² - 2*Vi*Vj*cos(θi - θj))
    totalLoss += G * (
      V[fromIdx] * V[fromIdx] + 
      V[toIdx] * V[toIdx] - 
      2 * V[fromIdx] * V[toIdx] * Math.cos(thetaIJ)
    );
  });
  
  return totalLoss;
}

// IMPROVED: Fitness function now includes losses in objective
function fitness(chromosome: number[], network: NetworkState, Ybus: Complex[][], vMin: number, vMax: number, alpha: number, beta: number, gamma: number): number {
  // Extract variables from chromosome
  let idx = 0;

  // Build voltage array (all buses)
  const V = new Array(network.nodes.length);
  network.nodes.forEach((node, i) => {
    if (node.type === NodeType.SLACK) {
      V[i] = node.voltage; // Fixed slack voltage
    } else {
      V[i] = chromosome[idx++]; // Optimized voltage
    }
  });

  // Build angle array (all buses)
  const theta = new Array(network.nodes.length);
  network.nodes.forEach((node, i) => {
    if (node.type === NodeType.SLACK) {
      theta[i] = 0; // Fixed slack angle
    } else if (node.type === NodeType.PV) {
      theta[i] = chromosome[idx++]; // Optimized angle for PV buses
    } else { // PQ buses
      theta[i] = chromosome[idx++]; // Optimized angle for PQ buses
    }
  });

  // Extract generation variables
  const slackPGens: number[] = [];
  const pvQGens: number[] = [];

  network.nodes.forEach(node => {
    if (node.type === NodeType.SLACK) {
      slackPGens.push(chromosome[idx++]);
    } else if (node.type === NodeType.PV) {
      pvQGens.push(chromosome[idx++]);
    }
  });

  // Calculate total system losses (objective to minimize)
  const losses = calculateTotalLosses(V, theta, network);

  let penalty = 0;
  let slackIdx = 0;
  let pvIdx = 0;

  network.nodes.forEach((bus, i) => {
    const calc = calculatePower(V, theta, Ybus, i);

    if (bus.type === NodeType.SLACK) {
      // Slack bus: P is optimized, Q is calculated
      const optimizedP = slackPGens[slackIdx++];
      const Pspec = optimizedP - bus.pLoad; // Use optimized P_gen
      penalty += alpha * Math.pow(calc.P - Pspec, 2);

      // Q constraint for slack (should be reasonable)
      const qViolation = Math.max(0, Math.abs(calc.Q) - 2.0);
      penalty += gamma * Math.pow(qViolation, 2);

    } else if (bus.type === NodeType.PV) {
      // PV bus: P is fixed, Q is optimized, V is fixed
      const Pspec = bus.pGen - bus.pLoad; // Fixed P
      penalty += alpha * Math.pow(calc.P - Pspec, 2);

      const optimizedQ = pvQGens[pvIdx++];
      const Qspec = optimizedQ - bus.qLoad; // Use optimized Q_gen
      penalty += beta * Math.pow(calc.Q - Qspec, 2);

      // Voltage should be close to target
      penalty += 1000 * Math.pow(V[i] - bus.voltage, 2);

    } else { // PQ bus
      // PQ bus: P and Q are fixed, V and θ are optimized
      const Pspec = bus.pGen - bus.pLoad;
      const Qspec = bus.qGen - bus.qLoad;
      penalty += alpha * Math.pow(calc.P - Pspec, 2);
      penalty += beta * Math.pow(calc.Q - Qspec, 2);
    }

    // Voltage constraint penalty (for all buses except slack)
    if (bus.type !== NodeType.SLACK && (V[i] < vMin || V[i] > vMax)) {
      penalty += 1000 * Math.pow(Math.max(vMin - V[i], V[i] - vMax, 0), 2);
    }
  });

  return losses + penalty;
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

// FIXED: Roulette wheel selection with proper fitness inversion
function rouletteWheel(population: number[][], fitnessScores: number[]): number[] {
  // For minimization: use inverse fitness
  const epsilon = 1e-10; // Avoid division by zero
  const invertedFitness = fitnessScores.map(f => 1 / (f + epsilon));
  
  const totalFitness = invertedFitness.reduce((sum, f) => sum + f, 0);
  const probabilities = invertedFitness.map(f => f / totalFitness);

  const cumulative = probabilities.reduce((acc, p, i) => {
    acc.push((acc[i - 1] || 0) + p);
    return acc;
  }, [] as number[]);

  const rand = Math.random();
  const selectedIndex = cumulative.findIndex(cum => cum >= rand);
  return population[selectedIndex !== -1 ? selectedIndex : population.length - 1];
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

// FIXED: Stochastic universal sampling with single random start
function stochasticUniversal(population: number[][], fitnessScores: number[]): number[][] {
  const numParents = 2;
  const epsilon = 1e-10;

  // Invert fitness for minimization
  const invertedFitness = fitnessScores.map(f => 1 / (f + epsilon));
  const totalFitness = invertedFitness.reduce((sum, f) => sum + f, 0);

  const probabilities = invertedFitness.map(f => f / totalFitness);
  const cumulative = probabilities.reduce((acc, p, i) => {
    acc.push((acc[i - 1] || 0) + p);
    return acc;
  }, [] as number[]);

  const selectedParents: number[][] = [];
  const spacing = 1.0 / numParents;
  
  // FIXED: Single random start point
  const start = Math.random() * spacing;

  for (let i = 0; i < numParents; i++) {
    const pointer = start + i * spacing;
    const selectedIndex = cumulative.findIndex(cum => cum >= pointer);
    const actualIndex = selectedIndex !== -1 ? selectedIndex : cumulative.length - 1;
    selectedParents.push([...population[actualIndex]]);
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

// ADDED: Gaussian random number generator
function gaussianRandom(): number {
  // Box-Muller transform for Gaussian distribution
  const u1 = Math.random();
  const u2 = Math.random();
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}

// IMPROVED: Mutation with Gaussian distribution for new chromosome structure
function mutate(chromosome: number[], rate: number, vMin: number, vMax: number, network: NetworkState): number[] {
  const mutated = [...chromosome];
  let idx = 0;

  const voltageSigma = 0.01; // Standard deviation for voltage mutation
  const angleSigma = 0.05;   // Standard deviation for angle mutation
  const pGenSigma = 1.0;     // Standard deviation for P generation mutation
  const qGenSigma = 0.5;     // Standard deviation for Q generation mutation

  // Mutate voltage magnitudes (PQ and PV buses)
  network.nodes.forEach(node => {
    if (node.type !== NodeType.SLACK) {
      if (Math.random() < rate) {
        mutated[idx] += gaussianRandom() * voltageSigma;
        mutated[idx] = Math.max(vMin, Math.min(vMax, mutated[idx]));
      }
      idx++;
    }
  });

  // Mutate voltage angles (PQ and PV buses)
  network.nodes.forEach(node => {
    if (node.type === NodeType.PQ || node.type === NodeType.PV) {
      if (Math.random() < rate) {
        mutated[idx] += gaussianRandom() * angleSigma;
      }
      idx++;
    }
  });

  // Mutate slack P generation
  network.nodes.forEach(node => {
    if (node.type === NodeType.SLACK) {
      if (Math.random() < rate) {
        mutated[idx] += gaussianRandom() * pGenSigma;
        // Keep within reasonable bounds (50% to 150% of base value)
        mutated[idx] = Math.max(node.pGen * 0.5, Math.min(node.pGen * 1.5, mutated[idx]));
      }
      idx++;
    }
  });

  // Mutate PV Q generation
  network.nodes.forEach(node => {
    if (node.type === NodeType.PV) {
      if (Math.random() < rate) {
        mutated[idx] += gaussianRandom() * qGenSigma;
        // Keep within reasonable bounds (±2.0 p.u.)
        mutated[idx] = Math.max(-2.0, Math.min(2.0, mutated[idx]));
      }
      idx++;
    }
  });

  return mutated;
}

// Simplified fitness function based on MATLAB approach
function fitnessMATLAB(X: number[], V1: number, Ybus: Complex[][], Pspec: number[], Qspec: number[],
                       Qg_min: number[], Qg_max: number[], alpha: number, beta: number, gamma: number,
                       network: NetworkState): { fitness: number, info: any } {
  const nb = network.nodes.length;
  // X = [V2 V3 ... V14] (voltages for buses 2-14)
  const V = new Array(nb);
  V[0] = V1; // Slack bus voltage fixed
  for (let i = 1; i < nb; i++) {
    V[i] = X[i-1]; // Optimized voltages
  }

  // Solve angles using Newton-Raphson (theta1 = 0 fixed)
  const theta0 = new Array(nb).fill(0);
  const theta = solveAnglesNewtonP(V, theta0, Ybus, Pspec);

  // Calculate P,Q from (V,theta)
  const { Pcalc, Qcalc } = calcPQ_MATLAB(V, theta, Ybus);

  // Active losses
  const Ploss = calcLosses_MATLAB(V, theta, network);

  // Mismatches
  const dP = Pspec.map((spec, i) => spec - Pcalc[i]);
  dP[0] = 0; // slack not penalized

  const dQ = Qspec.map((spec, i) => spec - Qcalc[i]);
  dQ[0] = 0; // slack not penalized
  // PV buses not penalized for Q
  network.nodes.forEach((node, i) => {
    if (node.type === 'PV') dQ[i] = 0;
  });

  // Qg violations for PV buses
  let penQg = 0;
  network.nodes.forEach((node, i) => {
    if (node.type === 'PV') {
      const Qg = Qcalc[i] + node.qLoad; // Qcalc + Qd
      if (Qg < Qg_min[i]) {
        penQg += Math.pow(Qg_min[i] - Qg, 2);
      } else if (Qg > Qg_max[i]) {
        penQg += Math.pow(Qg - Qg_max[i], 2);
      }
    }
  });

  const fitness = Ploss + alpha * dP.reduce((sum, dp) => sum + dp * dp, 0) +
                  beta * dQ.reduce((sum, dq) => sum + dq * dq, 0) + gamma * penQg;

  return {
    fitness,
    info: {
      Ploss,
      theta,
      Pcalc,
      Qcalc,
      dP,
      dQ
    }
  };
}

// Newton-Raphson to solve angles given voltages
function solveAnglesNewtonP(V: number[], theta0: number[], Ybus: Complex[][], Pspec: number[]): number[] {
  const nb = V.length;
  let theta = [...theta0];
  theta[0] = 0; // theta1 = 0 fixed

  // Initialize non-slack angles with small random values instead of zeros
  for (let i = 1; i < nb; i++) {
    theta[i] = (Math.random() - 0.5) * 0.1; // ±0.05 radians initial guess
  }

  const maxIt = 50; // Increased iterations
  const tol = 1e-6; // Relaxed tolerance
  const pqAngles = Array.from({length: nb-1}, (_, i) => i+1); // indices 1 to nb-1

  for (let it = 0; it < maxIt; it++) {
    const { Pcalc } = calcPQ_MATLAB(V, theta, Ybus);
    const mismatch = pqAngles.map(i => Pspec[i] - Pcalc[i]);

    if (Math.max(...mismatch.map(Math.abs)) < tol) {
      break;
    }

    // Jacobian dP/dtheta
    const J = jacobianP_theta_MATLAB(V, theta, Ybus, pqAngles);

    // Solve J * dth = mismatch
    const dth = solveLinearSystem(J, mismatch);
    if (!dth) break; // Singular matrix

    pqAngles.forEach((idx, i) => {
      theta[idx] = theta[idx] + dth[i];
    });
  }

  return theta;
}

// Jacobian for Newton method
function jacobianP_theta_MATLAB(V: number[], theta: number[], Ybus: Complex[][], idx: number[]): number[][] {
  const nb = V.length;
  const G = Ybus.map(row => row.map(y => y.real));
  const B = Ybus.map(row => row.map(y => y.imag));
  const m = idx.length;
  const J = Array(m).fill(0).map(() => Array(m).fill(0));

  // Calculate Q for diagonal terms
  const { Qcalc } = calcPQ_MATLAB(V, theta, Ybus);

  for (let a = 0; a < m; a++) {
    const i = idx[a];
    for (let b = 0; b < m; b++) {
      const k = idx[b];
      if (i === k) {
        // dP_i/dtheta_i = -Q_i - B_ii*V_i^2
        J[a][b] = -Qcalc[i] - B[i][i] * V[i] * V[i];
      } else {
        // dP_i/dtheta_k = V_i V_k (G_ik sin(θ_i-θ_k) - B_ik cos(θ_i-θ_k))
        const thik = theta[i] - theta[k];
        J[a][b] = V[i] * V[k] * (G[i][k] * Math.sin(thik) - B[i][k] * Math.cos(thik));
      }
    }
  }

  return J;
}

// Calculate P,Q injections
function calcPQ_MATLAB(V: number[], theta: number[], Ybus: Complex[][]): { Pcalc: number[], Qcalc: number[] } {
  const nb = V.length;
  const G = Ybus.map(row => row.map(y => y.real));
  const B = Ybus.map(row => row.map(y => y.imag));
  const Pcalc = new Array(nb).fill(0);
  const Qcalc = new Array(nb).fill(0);

  for (let i = 0; i < nb; i++) {
    for (let k = 0; k < nb; k++) {
      const th = theta[i] - theta[k];
      Pcalc[i] += V[i] * V[k] * (G[i][k] * Math.cos(th) + B[i][k] * Math.sin(th));
      Qcalc[i] += V[i] * V[k] * (G[i][k] * Math.sin(th) - B[i][k] * Math.cos(th));
    }
  }

  return { Pcalc, Qcalc };
}

// Calculate losses
function calcLosses_MATLAB(V: number[], theta: number[], network: NetworkState): number {
  let Ploss = 0;

  network.edges.forEach(edge => {
    const i = network.nodes.findIndex(node => node.id === edge.from);
    const j = network.nodes.findIndex(node => node.id === edge.to);
    if (i === -1 || j === -1) return;

    const R = edge.resistance;
    const X = edge.reactance;
    if (R === 0 && X === 0) return;

    const g = R / (R*R + X*X);
    const th = theta[i] - theta[j];
    Ploss += g * (V[i]*V[i] + V[j]*V[j] - 2*V[i]*V[j]*Math.cos(th));
  });

  return Ploss;
}

// Simple linear system solver (Gaussian elimination)
function solveLinearSystem(A: number[][], b: number[]): number[] | null {
  const n = A.length;
  const augmented = A.map((row, i) => [...row, b[i]]);

  // Forward elimination
  for (let i = 0; i < n; i++) {
    // Find pivot
    let maxRow = i;
    for (let k = i + 1; k < n; k++) {
      if (Math.abs(augmented[k][i]) > Math.abs(augmented[maxRow][i])) {
        maxRow = k;
      }
    }

    // Swap rows
    [augmented[i], augmented[maxRow]] = [augmented[maxRow], augmented[i]];

    // Check for singular matrix
    if (Math.abs(augmented[i][i]) < 1e-12) {
      return null;
    }

    // Eliminate
    for (let k = i + 1; k < n; k++) {
      const factor = augmented[k][i] / augmented[i][i];
      for (let j = i; j <= n; j++) {
        augmented[k][j] -= factor * augmented[i][j];
      }
    }
  }

  // Back substitution
  const x = new Array(n);
  for (let i = n - 1; i >= 0; i--) {
    x[i] = augmented[i][n];
    for (let j = i + 1; j < n; j++) {
      x[i] -= augmented[i][j] * x[j];
    }
    x[i] /= augmented[i][i];
  }

  return x;
}

// Main GA function - MATLAB inspired
export async function calculateLoadFlowGA(
  network: NetworkState,
  params: GAParameters,
  progressCallback?: (generation: number, fitness: number) => void
): Promise<GAResult> {
  const { populationSize, maxGenerations, crossoverRate, mutationRate, convergenceThreshold, minVoltage, maxVoltage } = params;

  // Build Y-bus matrix
  const Ybus = buildYbus(network);

  // Setup power specifications like MATLAB
  const nb = network.nodes.length;
  const Pspec = new Array(nb).fill(0);
  const Qspec = new Array(nb).fill(0);
  const Qg_min = new Array(nb).fill(0);
  const Qg_max = new Array(nb).fill(0);

  // Identify bus types
  const pv_buses: number[] = [];
  const pq_buses: number[] = [];

  network.nodes.forEach((node, i) => {
    if (node.type === 'PV') pv_buses.push(i);
    else if (node.type === 'PQ') pq_buses.push(i);
  });

  // Set specifications
  // For PV buses: Pspec = Pg - Pd (fixed)
  pv_buses.forEach(i => {
    Pspec[i] = network.nodes[i].pGen - network.nodes[i].pLoad;
  });

  // For PQ buses: Pspec = -Pd, Qspec = -Qd (loads are negative)
  pq_buses.forEach(i => {
    Pspec[i] = -network.nodes[i].pLoad;
    Qspec[i] = -network.nodes[i].qLoad;
  });

  // Qg limits for PV buses (simplified)
  pv_buses.forEach(i => {
    Qg_min[i] = -2.0; // Conservative limits
    Qg_max[i] = 2.0;
  });

  // GA parameters like MATLAB
  const D = nb - 1; // Variables: V2 to Vn (excluding slack voltage)

  // Initialize population
  let population: number[][] = [];
  for (let i = 0; i < populationSize; i++) {
    const chromosome = new Array(D);
    for (let j = 0; j < D; j++) {
      chromosome[j] = minVoltage + Math.random() * (maxVoltage - minVoltage);
    }
    population.push(chromosome);
  }

  const fitnessHistory: number[] = [];
  let bestSolution = null;
  let bestFitness = Infinity;
  let bestInfo: any = {};

  for (let gen = 0; gen < maxGenerations; gen++) {
    // Evaluate fitness
    const fitnessResults = population.map(ind =>
      fitnessMATLAB(ind, network.nodes[0].voltage, Ybus, Pspec, Qspec, Qg_min, Qg_max,
                    params.alpha, params.beta, params.gamma, network)
    );

    const fitnessScores = fitnessResults.map(r => r.fitness);
    const minFitness = Math.min(...fitnessScores);
    const minIdx = fitnessScores.indexOf(minFitness);

    if (minFitness < bestFitness) {
      bestFitness = minFitness;
      bestSolution = [...population[minIdx]];
      bestInfo = fitnessResults[minIdx].info;
    }

    fitnessHistory.push(bestFitness);

    if (progressCallback) {
      progressCallback(gen + 1, bestFitness);
    }

    // Check convergence
    if (bestFitness < convergenceThreshold) break;

    // New population (MATLAB style)
    const newPop: number[][] = [];
    const elite = 1; // Elitism

    // Elitism
    if (elite) {
      newPop.push([...bestSolution!]);
    }

    // Fill rest of population
    while (newPop.length < populationSize) {
      // Tournament selection (size 2)
      const idx1 = tournamentSelection(fitnessScores);
      const idx2 = tournamentSelection(fitnessScores);

      let p1 = [...population[idx1]];
      let p2 = [...population[idx2]];

      // Arithmetic crossover (like MATLAB)
      if (Math.random() < crossoverRate) {
        const lam = Math.random();
        const c1 = p1.map((v, i) => lam * v + (1 - lam) * p2[i]);
        const c2 = p1.map((v, i) => (1 - lam) * v + lam * p2[i]);
        p1 = c1;
        p2 = c2;
      }

      // Gaussian mutation
      p1 = gaussianMutation_MATLAB(p1, mutationRate, 0.01);
      p2 = gaussianMutation_MATLAB(p2, mutationRate, 0.01);

      // Enforce bounds
      p1 = p1.map(v => Math.min(Math.max(v, minVoltage), maxVoltage));
      p2 = p2.map(v => Math.min(Math.max(v, minVoltage), maxVoltage));

      newPop.push(p1);
      if (newPop.length < populationSize) {
        newPop.push(p2);
      }
    }

    population = newPop;

    // Allow UI to update
    if (gen % 10 === 0) {
      await new Promise(resolve => setTimeout(resolve, 0));
    }
  }

  // Reconstruct solution
  const fullVoltages = new Array(nb);
  fullVoltages[0] = network.nodes[0].voltage; // Slack voltage
  for (let i = 1; i < nb; i++) {
    fullVoltages[i] = bestSolution![i-1];
  }

  return {
    solution: {
      voltages: fullVoltages,
      angles: bestInfo.theta?.map((t: number) => t * 180 / Math.PI) || [],
      generations: {
        pGen: [], // Not optimized in this approach
        qGen: []  // Not optimized in this approach
      }
    },
    fitness: bestFitness,
    generations: fitnessHistory.length,
    fitnessHistory,
    converged: bestFitness < convergenceThreshold
  };
}

// Tournament selection (size 2)
function tournamentSelection(fitnessScores: number[]): number {
  const n = fitnessScores.length;
  const a = Math.floor(Math.random() * n);
  const b = Math.floor(Math.random() * n);
  return fitnessScores[a] < fitnessScores[b] ? a : b;
}

// Gaussian mutation (MATLAB style)
function gaussianMutation_MATLAB(chromosome: number[], Pm: number, sigma: number): number[] {
  return chromosome.map(v => {
    if (Math.random() < Pm) {
      return v + gaussianRandom() * sigma;
    }
    return v;
  });
}
