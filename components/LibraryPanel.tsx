
import React from 'react';
import { Zap, Activity, Home, Factory, Layers, Share2, Square } from 'lucide-react';
import { NodeType, ComponentCategory } from '../types';

export interface LibraryItem {
  id: string;
  name: string;
  type: NodeType;
  category: ComponentCategory;
  icon: React.ReactNode;
  defaultData: any;
}

interface LibraryPanelProps {
  onItemClick: (item: LibraryItem) => void;
}

const libraryItems: LibraryItem[] = [
  {
    id: 'bus-bar',
    name: 'Bus Bar',
    type: NodeType.PQ,
    category: 'Transmission',
    icon: <Square className="text-slate-800" size={18} fill="currentColor" fillOpacity={0.1} />,
    defaultData: { voltage: 1.0, pLoad: 0, qLoad: 0, name: 'New Bus Bar' }
  },
  {
    id: 'alternator-slack',
    name: 'Alternator (Slack)',
    type: NodeType.SLACK,
    category: 'Generation',
    icon: <Zap className="text-yellow-600" size={18} />,
    defaultData: { voltage: 1.05, pGen: 100, qGen: 0, name: 'Main Alternator' }
  },
  {
    id: 'alternator-pv',
    name: 'Alternator (PV)',
    type: NodeType.PV,
    category: 'Generation',
    icon: <Activity className="text-blue-600" size={18} />,
    defaultData: { voltage: 1.02, pGen: 50, qGen: 0, name: 'PV Alternator' }
  },
  {
    id: 'load-standard',
    name: 'Standard Load',
    type: NodeType.PQ,
    category: 'Load',
    icon: <Home className="text-emerald-600" size={18} />,
    defaultData: { pLoad: 40, qLoad: 15, name: 'Regional Load' }
  },
  {
    id: 'load-industrial',
    name: 'Industrial Load',
    type: NodeType.PQ,
    category: 'Load',
    icon: <Factory className="text-slate-700" size={18} />,
    defaultData: { pLoad: 150, qLoad: 60, name: 'Industrial Load' }
  },
  {
    id: 'branch-line',
    name: 'Transmission Line',
    type: NodeType.PQ,
    category: 'Transmission',
    icon: <Share2 className="text-slate-500" size={18} />,
    defaultData: { isTransformer: false, resistance: 0.02, reactance: 0.08, susceptance: 0.04, limit: 200, name: 'TX Line' }
  },
  {
    id: 'branch-transformer',
    name: 'Power Transformer',
    type: NodeType.PQ,
    category: 'Transmission',
    icon: <Layers className="text-orange-500" size={18} />,
    defaultData: { isTransformer: true, resistance: 0.005, reactance: 0.12, tapRatio: 1.02, limit: 150, name: 'Power Trafo' }
  }
];

const LibraryPanel: React.FC<LibraryPanelProps> = ({ onItemClick }) => {
  return (
    <div className="w-64 h-full bg-white border-r border-slate-200 flex flex-col overflow-hidden shadow-sm">
      <div className="p-4 border-b border-slate-100">
        <h2 className="font-bold text-slate-800 flex items-center gap-2">
          <Share2 size={18} className="text-blue-600" /> System Components
        </h2>
        <p className="text-[10px] text-slate-400 font-medium uppercase mt-1">Click to add to canvas</p>
      </div>
      
      <div className="flex-1 overflow-y-auto p-3 space-y-4">
        {['Generation', 'Load', 'Transmission'].map(category => (
          <div key={category} className="space-y-1">
            <h3 className="text-[10px] font-bold text-slate-400 uppercase px-2 mb-2 tracking-wider">{category}</h3>
            <div className="grid grid-cols-1 gap-1">
              {libraryItems.filter(i => i.category === category).map(item => (
                <button
                  key={item.id}
                  onClick={() => onItemClick(item)}
                  className="w-full flex items-center gap-3 p-2.5 rounded-lg border border-transparent hover:border-blue-200 hover:bg-blue-50 transition-all group text-left outline-none focus:ring-2 focus:ring-blue-100"
                >
                  <div className="p-2 bg-white rounded-md border border-slate-100 shadow-sm group-hover:shadow group-active:scale-95 transition-all">
                    {item.icon}
                  </div>
                  <span className="text-xs font-semibold text-slate-700">{item.name}</span>
                </button>
              ))}
            </div>
          </div>
        ))}
      </div>
      
      <div className="p-4 bg-slate-50 border-t border-slate-100">
        <div className="flex items-center gap-2 text-[10px] text-slate-500 font-medium italic leading-tight">
          <Activity size={12} className="shrink-0" />
          <span>New items appear in the center of the drawing area.</span>
        </div>
      </div>
    </div>
  );
};

export default LibraryPanel;
