
import { NetworkState, NodeType } from '../types';

export const IEEE_3_BUS: NetworkState = {
  nodes: [
    { id: '3b-1', name: 'Bus 1 (Slack)', type: NodeType.SLACK, category: 'Generation', voltage: 1.06, angle: 0, pGen: 259.49, qGen: -16.55, pLoad: 0, qLoad: 0, x: 500, y: 100, iconType: 'alternator-slack' },
    { id: '3b-2', name: 'Bus 2 (PV)', type: NodeType.PV, category: 'Generation', voltage: 1.045, angle: 0, pGen: 40, qGen: 0, pLoad: 21.7, qLoad: 12.7, x: 300, y: 200, iconType: 'alternator-pv' },
    { id: '3b-3', name: 'Bus 3 (PQ)', type: NodeType.PQ, category: 'Load', voltage: 1.0, angle: 0, pGen: 0, qGen: 0, pLoad: 94.2, qLoad: 19, x: 700, y: 300, iconType: 'load-standard' }
  ],
  edges: [
    { id: '3e-1', from: '3b-1', to: '3b-2', resistance: 0.02, reactance: 0.06, susceptance: 0.03, limit: 150, isTransformer: false },
    { id: '3e-2', from: '3b-2', to: '3b-3', resistance: 0.03, reactance: 0.09, susceptance: 0.04, limit: 150, isTransformer: false },
    { id: '3e-3', from: '3b-3', to: '3b-1', resistance: 0.02, reactance: 0.05, susceptance: 0.02, limit: 150, isTransformer: false }
  ]
};

export const IEEE_5_BUS: NetworkState = {
  nodes: [
    { id: '5b-1', name: 'Bus 1', type: NodeType.SLACK, category: 'Generation', voltage: 1.07, angle: 0, pGen: 0, qGen: 0, pLoad: 0, qLoad: 0, x: 400, y: 100, iconType: 'alternator-slack' },
    { id: '5b-2', name: 'Bus 2', type: NodeType.PV, category: 'Generation', voltage: 1.0, angle: 0, pGen: 40, qGen: 0, pLoad: 0, qLoad: 0, x: 150, y: 250, iconType: 'alternator-pv' },
    { id: '5b-3', name: 'Bus 3', type: NodeType.PQ, category: 'Load', voltage: 1.0, angle: 0, pGen: 0, qGen: 0, pLoad: 45, qLoad: 15, x: 250, y: 500, iconType: 'load-standard' },
    { id: '5b-4', name: 'Bus 4', type: NodeType.PQ, category: 'Load', voltage: 1.0, angle: 0, pGen: 0, qGen: 0, pLoad: 40, qLoad: 5, x: 550, y: 500, iconType: 'load-industrial' },
    { id: '5b-5', name: 'Bus 5', type: NodeType.PQ, category: 'Load', voltage: 1.0, angle: 0, pGen: 0, qGen: 0, pLoad: 60, qLoad: 10, x: 650, y: 250, iconType: 'load-industrial' }
  ],
  edges: [
    { id: '5e-1', from: '5b-1', to: '5b-2', resistance: 0.00281, reactance: 0.0281, susceptance: 0.00712, limit: 400 },
    { id: '5e-2', from: '5b-1', to: '5b-4', resistance: 0.00304, reactance: 0.0304, susceptance: 0.00658, limit: 426 },
    { id: '5e-3', from: '5b-1', to: '5b-5', resistance: 0.00064, reactance: 0.0064, susceptance: 0.03126, limit: 426 },
    { id: '5e-4', from: '5b-2', to: '5b-3', resistance: 0.00108, reactance: 0.0108, susceptance: 0.01852, limit: 426 },
    { id: '5e-5', from: '5b-3', to: '5b-4', resistance: 0.00297, reactance: 0.0297, susceptance: 0.00674, limit: 426 },
    { id: '5e-6', from: '5b-4', to: '5b-5', resistance: 0.00297, reactance: 0.0297, susceptance: 0.00674, limit: 240 }
  ]
};

export const IEEE_14_BUS: NetworkState = {
  nodes: [
    // Slack Bus (Generator)
    { id: '14b-1', name: 'Bus 1 (Slack)', type: NodeType.SLACK, category: 'Generation', voltage: 1.06, angle: 0, pGen: 0, qGen: 0, pLoad: 0, qLoad: 0, x: 500, y: 100, iconType: 'alternator-slack' },

    // PV Buses (Generators)
    { id: '14b-2', name: 'Bus 2 (PV)', type: NodeType.PV, category: 'Generation', voltage: 1.045, angle: 0, pGen: 40, qGen: 0, pLoad: 21.7, qLoad: 12.7, x: 300, y: 200, iconType: 'alternator-pv' },
    { id: '14b-3', name: 'Bus 3 (PV)', type: NodeType.PV, category: 'Generation', voltage: 1.01, angle: 0, pGen: 0, qGen: 0, pLoad: 94.2, qLoad: 19, x: 700, y: 200, iconType: 'alternator-pv' },
    { id: '14b-6', name: 'Bus 6 (PV)', type: NodeType.PV, category: 'Generation', voltage: 1.07, angle: 0, pGen: 0, qGen: 0, pLoad: 11.2, qLoad: 7.5, x: 800, y: 400, iconType: 'alternator-pv' },
    { id: '14b-8', name: 'Bus 8 (PV)', type: NodeType.PV, category: 'Generation', voltage: 1.09, angle: 0, pGen: 0, qGen: 0, pLoad: 0, qLoad: 0, x: 200, y: 500, iconType: 'alternator-pv' },

    // PQ Buses (Loads)
    { id: '14b-4', name: 'Bus 4 (PQ)', type: NodeType.PQ, category: 'Load', voltage: 1.0, angle: 0, pGen: 0, qGen: 0, pLoad: 47.8, qLoad: -3.9, x: 900, y: 300, iconType: 'load-standard' },
    { id: '14b-5', name: 'Bus 5 (PQ)', type: NodeType.PQ, category: 'Load', voltage: 1.0, angle: 0, pGen: 0, qGen: 0, pLoad: 7.6, qLoad: 1.6, x: 100, y: 300, iconType: 'load-standard' },
    { id: '14b-7', name: 'Bus 7 (PQ)', type: NodeType.PQ, category: 'Load', voltage: 1.0, angle: 0, pGen: 0, qGen: 0, pLoad: 0, qLoad: 0, x: 600, y: 500, iconType: 'bus-bar' },
    { id: '14b-9', name: 'Bus 9 (PQ)', type: NodeType.PQ, category: 'Load', voltage: 1.0, angle: 0, pGen: 0, qGen: 0, pLoad: 29.5, qLoad: 16.6, x: 400, y: 600, iconType: 'load-standard' },
    { id: '14b-10', name: 'Bus 10 (PQ)', type: NodeType.PQ, category: 'Load', voltage: 1.0, angle: 0, pGen: 0, qGen: 0, pLoad: 9, qLoad: 5.8, x: 600, y: 600, iconType: 'load-standard' },
    { id: '14b-11', name: 'Bus 11 (PQ)', type: NodeType.PQ, category: 'Load', voltage: 1.0, angle: 0, pGen: 0, qGen: 0, pLoad: 3.5, qLoad: 1.8, x: 200, y: 400, iconType: 'load-standard' },
    { id: '14b-12', name: 'Bus 12 (PQ)', type: NodeType.PQ, category: 'Load', voltage: 1.0, angle: 0, pGen: 0, qGen: 0, pLoad: 6.1, qLoad: 1.6, x: 800, y: 500, iconType: 'load-standard' },
    { id: '14b-13', name: 'Bus 13 (PQ)', type: NodeType.PQ, category: 'Load', voltage: 1.0, angle: 0, pGen: 0, qGen: 0, pLoad: 13.5, qLoad: 5.8, x: 100, y: 500, iconType: 'load-standard' },
    { id: '14b-14', name: 'Bus 14 (PQ)', type: NodeType.PQ, category: 'Load', voltage: 1.0, angle: 0, pGen: 0, qGen: 0, pLoad: 14.9, qLoad: 5, x: 900, y: 500, iconType: 'load-standard' }
  ],
  edges: [
    // Transmission lines
    { id: '14e-1', from: '14b-1', to: '14b-2', resistance: 0.01938, reactance: 0.05917, susceptance: 0.0528, limit: 100 },
    { id: '14e-2', from: '14b-1', to: '14b-5', resistance: 0.05403, reactance: 0.22304, susceptance: 0.0492, limit: 100 },
    { id: '14e-3', from: '14b-2', to: '14b-3', resistance: 0.04699, reactance: 0.19797, susceptance: 0.0438, limit: 100 },
    { id: '14e-4', from: '14b-2', to: '14b-4', resistance: 0.05811, reactance: 0.17632, susceptance: 0.034, limit: 100 },
    { id: '14e-5', from: '14b-2', to: '14b-5', resistance: 0.05695, reactance: 0.17388, susceptance: 0.0346, limit: 100 },
    { id: '14e-6', from: '14b-3', to: '14b-4', resistance: 0.06701, reactance: 0.17103, susceptance: 0.0128, limit: 100 },
    { id: '14e-7', from: '14b-4', to: '14b-5', resistance: 0.01335, reactance: 0.04211, susceptance: 0, limit: 100 },
    { id: '14e-8', from: '14b-4', to: '14b-7', resistance: 0, reactance: 0.20912, susceptance: 0, limit: 100, isTransformer: true },
    { id: '14e-9', from: '14b-4', to: '14b-9', resistance: 0, reactance: 0.55618, susceptance: 0, limit: 100, isTransformer: true },
    { id: '14e-10', from: '14b-5', to: '14b-6', resistance: 0, reactance: 0.25202, susceptance: 0, limit: 100, isTransformer: true },
    { id: '14e-11', from: '14b-6', to: '14b-11', resistance: 0.09498, reactance: 0.1989, susceptance: 0, limit: 100 },
    { id: '14e-12', from: '14b-6', to: '14b-12', resistance: 0.12291, reactance: 0.25581, susceptance: 0, limit: 100 },
    { id: '14e-13', from: '14b-6', to: '14b-13', resistance: 0.06615, reactance: 0.13027, susceptance: 0, limit: 100 },
    { id: '14e-14', from: '14b-7', to: '14b-8', resistance: 0, reactance: 0.17615, susceptance: 0, limit: 100, isTransformer: true },
    { id: '14e-15', from: '14b-7', to: '14b-9', resistance: 0.03181, reactance: 0.0845, susceptance: 0, limit: 100 },
    { id: '14e-16', from: '14b-9', to: '14b-10', resistance: 0.12711, reactance: 0.27038, susceptance: 0, limit: 100 },
    { id: '14e-17', from: '14b-9', to: '14b-14', resistance: 0.08205, reactance: 0.19207, susceptance: 0, limit: 100 },
    { id: '14e-18', from: '14b-10', to: '14b-11', resistance: 0.22092, reactance: 0.19988, susceptance: 0, limit: 100 },
    { id: '14e-19', from: '14b-12', to: '14b-13', resistance: 0.17093, reactance: 0.34802, susceptance: 0, limit: 100 },
    { id: '14e-20', from: '14b-13', to: '14b-14', resistance: 0.20746, reactance: 0.352, susceptance: 0, limit: 100 }
  ]
};

export const RADIAL_FEEDER: NetworkState = {
  nodes: [
    { id: 'r-1', name: 'Substation', type: NodeType.SLACK, category: 'Generation', voltage: 1.05, angle: 0, pGen: 10, qGen: 5, pLoad: 0, qLoad: 0, x: 400, y: 50, iconType: 'alternator-slack' },
    { id: 'r-2', name: 'Bus A', type: NodeType.PQ, category: 'Transmission', voltage: 1.0, angle: 0, pGen: 0, qGen: 0, pLoad: 2, qLoad: 1, x: 400, y: 200, iconType: 'bus-bar' },
    { id: 'r-3', name: 'Load B', type: NodeType.PQ, category: 'Load', voltage: 1.0, angle: 0, pGen: 0, qGen: 0, pLoad: 4, qLoad: 2, x: 250, y: 350, iconType: 'load-standard' },
    { id: 'r-4', name: 'Load C', type: NodeType.PQ, category: 'Load', voltage: 1.0, angle: 0, pGen: 0, qGen: 0, pLoad: 3, qLoad: 1.5, x: 550, y: 350, iconType: 'load-standard' }
  ],
  edges: [
    { id: 're-1', from: 'r-1', to: 'r-2', resistance: 0.01, reactance: 0.05, susceptance: 0, limit: 50, isTransformer: true, tapRatio: 1.05 },
    { id: 're-2', from: 'r-2', to: 'r-3', resistance: 0.05, reactance: 0.15, susceptance: 0, limit: 20 },
    { id: 're-3', from: 'r-2', to: 'r-4', resistance: 0.05, reactance: 0.15, susceptance: 0, limit: 20 }
  ]
};
