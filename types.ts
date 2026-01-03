
export enum NodeType {
  SLACK = 'SLACK',
  PQ = 'PQ',
  PV = 'PV'
}

export type ComponentCategory = 'Generation' | 'Load' | 'Renewable' | 'Transmission';

export interface BusNode {
  id: string;
  name: string;
  type: NodeType;
  category?: ComponentCategory;
  voltage: number; // p.u.
  angle: number; // degrees
  pGen: number; // MW
  qGen: number; // MVAr
  pLoad: number; // MW
  qLoad: number; // MVAr
  x: number;
  y: number;
  iconType?: string;
}

export interface BranchLine {
  id: string;
  from: string;
  to: string;
  resistance: number; // p.u.
  reactance: number; // p.u.
  susceptance: number; // p.u.
  limit: number; // MVA
  isTransformer?: boolean;
  tapRatio?: number;
  activeFlow?: number;
  reactiveFlow?: number;
  losses?: number;
}

export interface NetworkState {
  nodes: BusNode[];
  edges: BranchLine[];
}

export enum SelectionMethod {
  TOURNAMENT = 'TOURNAMENT',
  ROULETTE_WHEEL = 'ROULETTE_WHEEL',
  RANK_BASED = 'RANK_BASED',
  STOCHASTIC_UNIVERSAL = 'STOCHASTIC_UNIVERSAL'
}

export interface GAParameters {
  populationSize: number;
  maxGenerations: number;
  crossoverRate: number;
  mutationRate: number;
  convergenceThreshold: number;
  minVoltage: number;
  maxVoltage: number;
  selectionMethod: SelectionMethod;
  tournamentSize?: number; // For tournament selection
  alpha: number; // Penalty for delta P violations
  beta: number; // Penalty for delta Q violations
  gamma: number; // Penalty for Qg violations
}

export interface GAResult {
  solution: {
    voltages: number[];
    angles: number[];
    generations?: {
      pGen: number[];  // Optimized active power generation for slack buses
      qGen: number[];  // Optimized reactive power generation for PV buses
    };
  };
  fitness: number;
  generations: number;
  fitnessHistory: number[];
  converged: boolean;
}

export interface CalculationResult {
  nodes: {
    id: string;
    voltage: number;
    angle: number;
    pGen: number;
    qGen: number;
  }[];
  edges: {
    id: string;
    pFlow: number;
    qFlow: number;
    losses: number;
  }[];
  summary: string;
  gaResult?: GAResult;
}
