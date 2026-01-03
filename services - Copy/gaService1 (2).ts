// ============= FIX 1: Y-bus Construction =============
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

// ============= FIX 2: Calculate Total Losses =============
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

// ============= FIX 3: Improved Fitness Function =============
function fitness(
  chromosome: number[], 
  network: NetworkState, 
  Ybus: Complex[][], 
  vMin: number, 
  vMax: number, 
  alpha: number, 
  beta: number, 
  gamma: number
): number {
  const n = network.nodes.length;
  const V = chromosome.slice(0, n);
  const theta = chromosome.slice(n);

  // Calculate total system losses (objective to minimize)
  const losses = calculateTotalLosses(V, theta, network);
  
  let penalty = 0;

  network.nodes.forEach((bus, i) => {
    if (bus.type === NodeType.SLACK) return;

    const calc = calculatePower(V, theta, Ybus, i);
    const Pspec = bus.pGen - bus.pLoad;
    const Qspec = bus.qGen - bus.qLoad;

    // Power balance penalties
    penalty += alpha * Math.pow(calc.P - Pspec, 2);
    
    if (bus.type === NodeType.PQ) {
      penalty += beta * Math.pow(calc.Q - Qspec, 2);
    }

    // Voltage constraint penalty
    if (V[i] < vMin || V[i] > vMax) {
      penalty += 1000 * Math.pow(Math.max(vMin - V[i], V[i] - vMax, 0), 2);
    }

    // FIXED: Reactive power limits with actual constraints
    if (bus.type === NodeType.PV || bus.type === NodeType.SLACK) {
      const Qg = calc.Q + bus.qLoad; // Qg = Qcalc + Qload
      
      // Use actual limits from bus data (if available)
      const qMin = bus.qMin ?? -Infinity;
      const qMax = bus.qMax ?? Infinity;
      
      if (Qg < qMin) {
        penalty += gamma * Math.pow(qMin - Qg, 2);
      } else if (Qg > qMax) {
        penalty += gamma * Math.pow(Qg - qMax, 2);
      }
    }
  });

  return losses + penalty;
}

// ============= FIX 4: Roulette Wheel Selection =============
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

// ============= FIX 5: Stochastic Universal Sampling =============
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

// ============= FIX 6: Gaussian Mutation =============
function gaussianRandom(): number {
  // Box-Muller transform for Gaussian distribution
  const u1 = Math.random();
  const u2 = Math.random();
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}

function mutate(
  chromosome: number[], 
  rate: number, 
  vMin: number, 
  vMax: number, 
  network: NetworkState
): number[] {
  const n = network.nodes.length;
  const mutated = [...chromosome];
  
  const voltageSigma = 0.01; // Standard deviation for voltage mutation
  const angleSigma = 0.05;   // Standard deviation for angle mutation

  for (let i = 0; i < chromosome.length; i++) {
    if (Math.random() < rate) {
      if (i < n) {
        // Voltage mutation (Gaussian)
        if (network.nodes[i].type !== NodeType.SLACK) {
          mutated[i] += gaussianRandom() * voltageSigma;
          mutated[i] = Math.max(vMin, Math.min(vMax, mutated[i]));
        }
      } else {
        // Angle mutation (Gaussian)
        const busIdx = i - n;
        if (network.nodes[busIdx].type !== NodeType.SLACK) {
          mutated[i] += gaussianRandom() * angleSigma;
        }
      }
    }
  }

  return mutated;
}

// ============= CRITICAL NOTE =============
/*
  FUNDAMENTAL ISSUE: This GA optimizes (V, θ) directly without ensuring
  they satisfy power flow equations. For proper power flow:
  
  OPTION A: Only optimize V for PQ buses, solve θ using Newton-Raphson
  OPTION B: Run full Newton-Raphson inside fitness function
  
  Current approach will never converge to valid power flow solution
  because it treats power flow equations as soft constraints (penalties)
  rather than hard constraints.
  
  Recommendation: Implement Option A for better convergence.
*/