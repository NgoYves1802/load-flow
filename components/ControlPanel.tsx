
import React, { useState, useEffect } from 'react';
import { BusNode, BranchLine, NodeType, CalculationResult, GAParameters, GAResult, SelectionMethod } from '../types';
import { Settings, Play, Info, Trash2, Cpu, AlertTriangle, Activity, Layers, Share2, Zap, Home, Sliders, BarChart3 } from 'lucide-react';

interface Props {
  selectedId: string | null;
  selectedType: 'node' | 'edge' | null;
  nodes: BusNode[];
  edges: BranchLine[];
  onUpdateNode: (node: BusNode) => void;
  onUpdateEdge: (edge: BranchLine) => void;
  onDeleteNode: (id: string) => void;
  onDeleteEdge: (id: string) => void;
  onRunLoadFlow: () => void;
  isCalculating: boolean;
  results: CalculationResult | null;
  gaParams: GAParameters;
  onUpdateGAParams: (params: GAParameters) => void;
}

const ControlPanel: React.FC<Props> = ({
  selectedId, selectedType, nodes, edges, onUpdateNode, onUpdateEdge, onDeleteNode, onDeleteEdge, onRunLoadFlow, isCalculating, results, gaParams, onUpdateGAParams
}) => {
  const [activeTab, setActiveTab] = useState<'parameters' | 'characteristics'>('characteristics');

  const selectedNode = selectedType === 'node' ? nodes.find(n => n.id === selectedId) : null;
  const selectedEdge = selectedType === 'edge' ? edges.find(e => e.id === selectedId) : null;

  const handleNumChange = (val: string, min: number = -Infinity): number => {
    const parsed = parseFloat(val);
    if (isNaN(parsed)) return 0;
    return Math.max(min, parsed);
  };

  const getBusHeading = (node: BusNode) => {
    if (node.type === NodeType.SLACK || node.type === NodeType.PV) return "Alternator Settings";
    if (node.pLoad > 0 || node.qLoad > 0) return "Load Characteristics";
    return "Bus Bar Configuration";
  };

  return (
    <div className="w-80 h-full bg-white border-l border-slate-200 flex flex-col overflow-hidden">
      <div className="p-4 border-b border-slate-100 flex items-center justify-between">
        <h2 className="font-bold text-slate-800 flex items-center gap-2">
          <Settings size={18} /> Characteristics
        </h2>
        <button 
          onClick={onRunLoadFlow}
          disabled={isCalculating || nodes.length === 0}
          className={`flex items-center gap-2 px-4 py-2 rounded-lg font-bold text-white transition-all shadow-md
            ${isCalculating ? 'bg-slate-400 cursor-not-allowed' : 'bg-blue-600 hover:bg-blue-700 active:scale-95'}
          `}
        >
          {isCalculating ? <Cpu className="animate-spin" size={18} /> : <Play size={18} />}
          {isCalculating ? 'Solving...' : 'Run Solver'}
        </button>
      </div>

      {/* Tab Navigation */}
      <div className="flex border-b border-slate-200">
        <button
          onClick={() => setActiveTab('parameters')}
          className={`flex-1 py-3 px-2 text-xs font-medium transition-colors flex items-center justify-center gap-1 ${
            activeTab === 'parameters'
              ? 'text-blue-600 border-b-2 border-blue-600 bg-blue-50'
              : 'text-slate-500 hover:text-slate-700 hover:bg-slate-50'
          }`}
        >
          <Sliders size={14} />
          Parameters
        </button>
        <button
          onClick={() => setActiveTab('characteristics')}
          className={`flex-1 py-3 px-2 text-xs font-medium transition-colors flex items-center justify-center gap-1 ${
            activeTab === 'characteristics'
              ? 'text-blue-600 border-b-2 border-blue-600 bg-blue-50'
              : 'text-slate-500 hover:text-slate-700 hover:bg-slate-50'
          }`}
        >
          <Settings size={14} />
          Characteristics
        </button>
      </div>

      <div className="flex-1 overflow-y-auto p-4 space-y-6 scrollbar-hide">
        {activeTab === 'parameters' && (
          <div className="space-y-6">
            {/* GA Parameters */}
            <div className="space-y-4">
              <h3 className="text-sm font-bold text-slate-800 flex items-center gap-2">
                <Cpu size={16} />
                Genetic Algorithm Parameters
              </h3>

              <div className="grid grid-cols-2 gap-4">
                <div className="space-y-2">
                  <label className="text-xs font-medium text-slate-600">Population Size</label>
                  <input
                    type="number"
                    value={gaParams.populationSize}
                    onChange={e => onUpdateGAParams({...gaParams, populationSize: parseInt(e.target.value) || 100})}
                    min="20"
                    max="1000"
                    className="w-full px-3 py-2 border border-slate-200 rounded-md text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                  />
                </div>
                <div className="space-y-2">
                  <label className="text-xs font-medium text-slate-600">Max Generations</label>
                  <input
                    type="number"
                    value={gaParams.maxGenerations}
                    onChange={e => onUpdateGAParams({...gaParams, maxGenerations: parseInt(e.target.value) || 300})}
                    min="50"
                    max="2000"
                    className="w-full px-3 py-2 border border-slate-200 rounded-md text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                  />
                </div>
                <div className="space-y-2">
                  <label className="text-xs font-medium text-slate-600">Crossover Rate</label>
                  <input
                    type="number"
                    value={gaParams.crossoverRate}
                    onChange={e => onUpdateGAParams({...gaParams, crossoverRate: parseFloat(e.target.value) || 0.8})}
                    min="0"
                    max="1"
                    step="0.05"
                    className="w-full px-3 py-2 border border-slate-200 rounded-md text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                  />
                </div>
                <div className="space-y-2">
                  <label className="text-xs font-medium text-slate-600">Mutation Rate</label>
                  <input
                    type="number"
                    value={gaParams.mutationRate}
                    onChange={e => onUpdateGAParams({...gaParams, mutationRate: parseFloat(e.target.value) || 0.05})}
                    min="0"
                    max="1"
                    step="0.01"
                    className="w-full px-3 py-2 border border-slate-200 rounded-md text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                  />
                </div>
                <div className="space-y-2">
                  <label className="text-xs font-medium text-slate-600">Convergence Threshold</label>
                  <input
                    type="number"
                    value={gaParams.convergenceThreshold}
                    onChange={e => onUpdateGAParams({...gaParams, convergenceThreshold: parseFloat(e.target.value) || 0.001})}
                    min="0"
                    step="0.0001"
                    className="w-full px-3 py-2 border border-slate-200 rounded-md text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                  />
                </div>
                <div className="space-y-2">
                  <label className="text-xs font-medium text-slate-600">Selection Method</label>
                  <select
                    value={gaParams.selectionMethod}
                    onChange={e => onUpdateGAParams({...gaParams, selectionMethod: e.target.value as SelectionMethod})}
                    className="w-full px-3 py-2 border border-slate-200 rounded-md text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                  >
                    <option value={SelectionMethod.TOURNAMENT}>Tournament</option>
                    <option value={SelectionMethod.ROULETTE_WHEEL}>Roulette Wheel</option>
                    <option value={SelectionMethod.RANK_BASED}>Rank Based</option>
                    <option value={SelectionMethod.STOCHASTIC_UNIVERSAL}>Stochastic Universal</option>
                  </select>
                </div>
              </div>

              {/* Tournament Size - only show for tournament selection */}
              {gaParams.selectionMethod === SelectionMethod.TOURNAMENT && (
                <div className="space-y-2">
                  <label className="text-xs font-medium text-slate-600">Tournament Size</label>
                  <input
                    type="number"
                    value={gaParams.tournamentSize || 5}
                    onChange={e => onUpdateGAParams({...gaParams, tournamentSize: parseInt(e.target.value) || 5})}
                    min="2"
                    max="10"
                    className="w-full px-3 py-2 border border-slate-200 rounded-md text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                  />
                </div>
              )}
            </div>

            {/* Penalty Parameters */}
            <div className="space-y-4">
              <h3 className="text-sm font-bold text-slate-800 flex items-center gap-2">
                <AlertTriangle size={16} />
                Penalty Weights
              </h3>

              <div className="grid grid-cols-1 gap-4">
                <div className="space-y-2">
                  <label className="text-xs font-medium text-slate-600">Alpha (δP violations)</label>
                  <input
                    type="number"
                    value={gaParams.alpha}
                    onChange={e => onUpdateGAParams({...gaParams, alpha: parseFloat(e.target.value) || 1e6})}
                    min="0"
                    step="1000"
                    className="w-full px-3 py-2 border border-slate-200 rounded-md text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                  />
                  <p className="text-[10px] text-slate-500">Penalty for active power balance violations</p>
                </div>

                <div className="space-y-2">
                  <label className="text-xs font-medium text-slate-600">Beta (δQ violations)</label>
                  <input
                    type="number"
                    value={gaParams.beta}
                    onChange={e => onUpdateGAParams({...gaParams, beta: parseFloat(e.target.value) || 1e6})}
                    min="0"
                    step="1000"
                    className="w-full px-3 py-2 border border-slate-200 rounded-md text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                  />
                  <p className="text-[10px] text-slate-500">Penalty for reactive power balance violations</p>
                </div>

                <div className="space-y-2">
                  <label className="text-xs font-medium text-slate-600">Gamma (Qg violations)</label>
                  <input
                    type="number"
                    value={gaParams.gamma}
                    onChange={e => onUpdateGAParams({...gaParams, gamma: parseFloat(e.target.value) || 1e5})}
                    min="0"
                    step="1000"
                    className="w-full px-3 py-2 border border-slate-200 rounded-md text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                  />
                  <p className="text-[10px] text-slate-500">Penalty for generator reactive power limit violations</p>
                </div>
              </div>
            </div>

            {/* Voltage Constraints */}
            <div className="space-y-4">
              <h3 className="text-sm font-bold text-slate-800 flex items-center gap-2">
                <Zap size={16} />
                Voltage Constraints
              </h3>
              <div className="grid grid-cols-2 gap-4">
                <div className="space-y-2">
                  <label className="text-xs font-medium text-slate-600">Min Voltage (p.u.)</label>
                  <input
                    type="number"
                    value={gaParams.minVoltage}
                    onChange={e => onUpdateGAParams({...gaParams, minVoltage: parseFloat(e.target.value) || 0.95})}
                    min="0.8"
                    max="1.0"
                    step="0.01"
                    className="w-full px-3 py-2 border border-slate-200 rounded-md text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                  />
                </div>
                <div className="space-y-2">
                  <label className="text-xs font-medium text-slate-600">Max Voltage (p.u.)</label>
                  <input
                    type="number"
                    value={gaParams.maxVoltage}
                    onChange={e => onUpdateGAParams({...gaParams, maxVoltage: parseFloat(e.target.value) || 1.05})}
                    min="1.0"
                    max="1.2"
                    step="0.01"
                    className="w-full px-3 py-2 border border-slate-200 rounded-md text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                  />
                </div>
              </div>
            </div>
          </div>
        )}

        {activeTab === 'characteristics' && (
          <div className="space-y-6">
            {/* Component Editing */}
            {!selectedId && (
              <div className="flex flex-col items-center justify-center h-64 text-slate-400 text-center">
                <Info size={40} className="mb-4 opacity-20" />
                <p className="text-sm">Select a component on the canvas to edit its electrical characteristics.</p>
              </div>
            )}

            {selectedNode && (
              <div className="space-y-4">
                <div className="flex justify-between items-center">
                  <span className="text-xs font-bold uppercase tracking-wider text-slate-400">
                    {getBusHeading(selectedNode)}
                  </span>
                  <button onClick={() => onDeleteNode(selectedNode.id)} className="text-red-500 hover:text-red-700">
                    <Trash2 size={16} />
                  </button>
                </div>

                <div className="space-y-2">
                  <label className="text-xs font-medium text-slate-500">Tag / Label</label>
                  <input
                    value={selectedNode.name}
                    onChange={e => onUpdateNode({...selectedNode, name: e.target.value})}
                    className="w-full px-3 py-2 border border-slate-200 rounded-md text-sm focus:ring-2 focus:ring-blue-500 outline-none transition-all"
                  />
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-1">
                    <label className="text-[10px] font-bold text-slate-400 uppercase">Target Voltage (V)</label>
                    <input
                      type="number" step="0.001" min="0"
                      value={selectedNode.voltage}
                      onChange={e => onUpdateNode({...selectedNode, voltage: handleNumChange(e.target.value, 0)})}
                      className="w-full px-2 py-1 border border-slate-200 rounded text-sm fira-code focus:border-blue-400 outline-none"
                    />
                  </div>
                  <div className="space-y-1">
                    <label className="text-[10px] font-bold text-slate-400 uppercase">Phase Angle (δ)</label>
                    <input
                      type="number" step="0.1"
                      value={selectedNode.angle}
                      disabled={selectedNode.type === NodeType.SLACK}
                      onChange={e => onUpdateNode({...selectedNode, angle: handleNumChange(e.target.value)})}
                      className={`w-full px-2 py-1 border border-slate-200 rounded text-sm fira-code outline-none ${selectedNode.type === NodeType.SLACK ? 'bg-slate-50 text-slate-400' : 'focus:border-blue-400'}`}
                    />
                  </div>
                </div>

                {(selectedNode.type !== NodeType.PQ || selectedNode.pGen > 0) && (
                  <div className="space-y-3 pt-2 border-t border-slate-100">
                    <h3 className="text-[10px] font-bold text-slate-400 uppercase flex items-center gap-2">
                      <Zap size={12} /> Alternator Parameters
                    </h3>
                    <div className="grid grid-cols-2 gap-3">
                      <div className="space-y-1">
                        <label className="text-[10px] text-slate-500">Active Power (P)</label>
                        <input
                          type="number" step="0.1"
                          value={selectedNode.pGen}
                          onChange={e => onUpdateNode({...selectedNode, pGen: handleNumChange(e.target.value)})}
                          className="w-full px-2 py-1 border border-slate-200 rounded text-xs fira-code focus:border-blue-400 outline-none"
                        />
                      </div>
                      <div className="space-y-1">
                        <label className="text-[10px] text-slate-500">Reactive Power (Q)</label>
                        <input
                          type="number" step="0.1"
                          value={selectedNode.qGen}
                          onChange={e => onUpdateNode({...selectedNode, qGen: handleNumChange(e.target.value)})}
                          className="w-full px-2 py-1 border border-slate-200 rounded text-xs fira-code focus:border-blue-400 outline-none"
                        />
                      </div>
                    </div>
                  </div>
                )}

                {(selectedNode.type === NodeType.PQ || selectedNode.pLoad > 0) && (
                  <div className="space-y-3 pt-2 border-t border-slate-100">
                    <h3 className="text-[10px] font-bold text-slate-400 uppercase flex items-center gap-2">
                      <Home size={12} /> Load Demand
                    </h3>
                    <div className="grid grid-cols-2 gap-3">
                      <div className="space-y-1">
                        <label className="text-[10px] text-slate-500">Active Load (P)</label>
                        <input
                          type="number" step="0.1" min="0"
                          value={selectedNode.pLoad}
                          onChange={e => onUpdateNode({...selectedNode, pLoad: handleNumChange(e.target.value, 0)})}
                          className="w-full px-2 py-1 border border-slate-200 rounded text-xs fira-code focus:border-blue-400 outline-none"
                        />
                      </div>
                      <div className="space-y-1">
                        <label className="text-[10px] text-slate-500">Reactive Load (Q)</label>
                        <input
                          type="number" step="0.1" min="0"
                          value={selectedNode.qLoad}
                          onChange={e => onUpdateNode({...selectedNode, qLoad: handleNumChange(e.target.value, 0)})}
                          className="w-full px-2 py-1 border border-slate-200 rounded text-xs fira-code focus:border-blue-400 outline-none"
                        />
                      </div>
                    </div>
                  </div>
                )}
              </div>
            )}

            {selectedEdge && (
              <div className="space-y-4">
                <div className="flex justify-between items-center">
                  <span className="text-xs font-bold uppercase tracking-wider text-slate-400">
                    {selectedEdge.isTransformer ? "Power Transformer" : "Transmission Line"}
                  </span>
                  <button onClick={() => onDeleteEdge(selectedEdge.id)} className="text-red-500 hover:text-red-700">
                    <Trash2 size={16} />
                  </button>
                </div>

                <div className="flex items-center gap-2 p-3 bg-slate-50 rounded-lg border border-slate-100">
                  <input
                    type="checkbox"
                    checked={!!selectedEdge.isTransformer}
                    onChange={e => onUpdateEdge({...selectedEdge, isTransformer: e.target.checked})}
                    id="isTransformer"
                    className="w-4 h-4 text-blue-600 rounded"
                  />
                  <label htmlFor="isTransformer" className="text-xs font-semibold text-slate-700 flex items-center gap-2 cursor-pointer">
                    <Layers size={14} className="text-orange-500" /> Type: Transformer
                  </label>
                </div>

                {selectedEdge.isTransformer && (
                  <div className="space-y-1">
                    <label className="text-[10px] text-slate-500 uppercase font-bold">Tap Setting (p.u.)</label>
                    <input
                      type="number" step="0.01" min="0.5" max="1.5"
                      value={selectedEdge.tapRatio || 1.0}
                      onChange={e => onUpdateEdge({...selectedEdge, tapRatio: handleNumChange(e.target.value, 0.5)})}
                      className="w-full px-2 py-1 border border-slate-200 rounded text-sm fira-code focus:border-blue-400 outline-none"
                    />
                  </div>
                )}

                <div className="grid grid-cols-2 gap-4">
                  <div className="space-y-1">
                    <label className="text-[10px] text-slate-500">Resistance (R)</label>
                    <input
                      type="number" step="0.00001" min="0"
                      value={selectedEdge.resistance}
                      onChange={e => onUpdateEdge({...selectedEdge, resistance: handleNumChange(e.target.value, 0)})}
                      className="w-full px-2 py-1 border border-slate-200 rounded text-sm fira-code focus:border-blue-400 outline-none"
                    />
                  </div>
                  <div className="space-y-1">
                    <label className="text-[10px] text-slate-500">Reactance (X)</label>
                    <input
                      type="number" step="0.00001" min="0"
                      value={selectedEdge.reactance}
                      onChange={e => onUpdateEdge({...selectedEdge, reactance: handleNumChange(e.target.value, 0)})}
                      className="w-full px-2 py-1 border border-slate-200 rounded text-sm fira-code focus:border-blue-400 outline-none"
                    />
                  </div>
                </div>

                {!selectedEdge.isTransformer && (
                  <div className="space-y-1">
                    <label className="text-[10px] text-slate-500">Line Charging (B/2)</label>
                    <input
                      type="number" step="0.001" min="0"
                      value={selectedEdge.susceptance}
                      onChange={e => onUpdateEdge({...selectedEdge, susceptance: handleNumChange(e.target.value, 0)})}
                      className="w-full px-2 py-1 border border-slate-200 rounded text-sm fira-code focus:border-blue-400 outline-none"
                    />
                  </div>
                )}

                <div className="space-y-2">
                  <label className="text-xs font-medium text-slate-500">Thermal Limit (MVA)</label>
                  <input
                    type="number" step="1" min="0"
                    value={selectedEdge.limit}
                    onChange={e => onUpdateEdge({...selectedEdge, limit: handleNumChange(e.target.value, 0)})}
                    className="w-full px-3 py-2 border border-slate-200 rounded-md text-sm focus:ring-2 focus:ring-blue-500 outline-none transition-all"
                  />
                </div>
              </div>
            )}
          </div>
        )}

      </div>
    </div>
  );
};

export default ControlPanel;
