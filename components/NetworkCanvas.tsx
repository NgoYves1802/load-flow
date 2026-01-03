
import React, { useState, useRef } from 'react';
import { BusNode, BranchLine, NodeType } from '../types';
import { Zap, Activity, Home, Trash2, Plus, Square, Share2, Layers, Minus, RotateCcw } from 'lucide-react';

interface Props {
  nodes: BusNode[];
  edges: BranchLine[];
  onAddNode: (node: BusNode) => void;
  onAddEdge: (edge: BranchLine) => void;
  onUpdateNode: (node: BusNode, isTransient?: boolean) => void;
  onUpdateEdge: (edge: BranchLine) => void;
  onDeleteNode: (id: string) => void;
  onDeleteEdge: (id: string) => void;
  selectedId: string | null;
  onSelect: (id: string | null, type: 'node' | 'edge' | null) => void;
  zoom?: number;
  onZoomChange?: (zoom: number) => void;
}

const NetworkCanvas: React.FC<Props> = ({
  nodes, edges, onAddNode, onAddEdge, onUpdateNode, onUpdateEdge, onDeleteNode, onDeleteEdge, selectedId, onSelect, zoom = 1, onZoomChange
}) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const [isLinking, setIsLinking] = useState<string | null>(null);

  const handleZoomIn = () => {
    const newZoom = Math.min(zoom * 1.2, 5); // Max zoom 5x
    onZoomChange?.(newZoom);
  };

  const handleZoomOut = () => {
    const newZoom = Math.max(zoom / 1.2, 0.1); // Min zoom 0.1x
    onZoomChange?.(newZoom);
  };

  const handleResetZoom = () => {
    onZoomChange?.(1);
  };

  const handleCanvasClick = (e: React.MouseEvent) => {
    if (e.target === containerRef.current) {
      onSelect(null, null);
      setIsLinking(null);
    }
  };

  const handleNodeClick = (id: string, e: React.MouseEvent) => {
    e.stopPropagation();
    if (isLinking && isLinking !== id) {
      const newEdge: BranchLine = {
        id: `L-${Date.now()}`,
        from: isLinking,
        to: id,
        resistance: 0.05,
        reactance: 0.15,
        susceptance: 0.02,
        limit: 100,
        isTransformer: false
      };
      onAddEdge(newEdge);
      setIsLinking(null);
      onSelect(newEdge.id, 'edge');
    } else {
      onSelect(id, 'node');
    }
  };

  const getNodeIcon = (node: BusNode) => {
    if (node.iconType?.includes('alternator')) return <Zap className="text-yellow-500" size={20} />;
    if (node.iconType?.includes('load')) return <Home className="text-emerald-500" size={20} />;
    if (node.iconType === 'bus-bar') return <Square className="text-slate-700" size={20} fill="currentColor" fillOpacity={0.1} />;
    
    switch (node.type) {
      case NodeType.SLACK: return <Zap className="text-yellow-500" size={20} />;
      case NodeType.PV: return <Activity className="text-blue-500" size={20} />;
      case NodeType.PQ: return <Square className="text-slate-400" size={20} />;
      default: return <Activity size={20} />;
    }
  };

  return (
    <div className="relative w-full h-full bg-white overflow-auto rounded-xl border border-slate-200 transition-colors duration-200 shadow-inner" data-canvas="network-canvas">
      <div
        ref={containerRef}
        className="relative cursor-crosshair bg-gradient-to-br from-slate-50 to-white"
        onClick={handleCanvasClick}
        style={{
          width: '4000px',
          height: '3000px',
          transform: `scale(${zoom})`,
          transformOrigin: 'top left',
          backgroundImage: `
            linear-gradient(rgba(0,0,0,0.1) 1px, transparent 1px),
            linear-gradient(90deg, rgba(0,0,0,0.1) 1px, transparent 1px)
          `,
          backgroundSize: `${20 * zoom}px ${20 * zoom}px`,
        }}
        data-zoom-level={zoom}
      >
        <svg
          className="absolute inset-0 pointer-events-none"
          width="4000"
          height="3000"
          style={{ overflow: 'visible' }}
        >
          {edges.map(edge => {
            const fromNode = nodes.find(n => n.id === edge.from);
            const toNode = nodes.find(n => n.id === edge.to);
            if (!fromNode || !toNode) return null;
            
            const isSelected = selectedId === edge.id;
            const midX = (fromNode.x + toNode.x) / 2;
            const midY = (fromNode.y + toNode.y) / 2;
            
            return (
              <g key={edge.id} className="pointer-events-auto cursor-pointer" onClick={() => onSelect(edge.id, 'edge')}>
                <line 
                  x1={fromNode.x} y1={fromNode.y} 
                  x2={toNode.x} y2={toNode.y} 
                  stroke={isSelected ? '#3b82f6' : edge.isTransformer ? '#f59e0b' : '#94a3b8'} 
                  strokeWidth={isSelected ? 4 : edge.isTransformer ? 3 : 2}
                  strokeDasharray={edge.isTransformer ? "6,4" : "0"}
                />
                {edge.isTransformer && (
                  <circle cx={midX} cy={midY} r="8" fill="white" stroke="#f59e0b" strokeWidth="2" />
                )}
                <text 
                  x={midX} 
                  y={midY - 15}
                  className="text-[10px] fill-slate-500 font-bold pointer-events-none select-none"
                  textAnchor="middle"
                >
                  {edge.activeFlow ? `${edge.activeFlow.toFixed(1)} MW` : ''}
                </text>
              </g>
            );
          })}
        </svg>

        {nodes.map(node => (
          <div 
            key={node.id}
            onClick={(e) => handleNodeClick(node.id, e)}
            onMouseDown={(e) => {
              if (e.button !== 0) return;
              const startX = e.clientX;
              const startY = e.clientY;
              const initialNodeX = node.x;
              const initialNodeY = node.y;
              let currentX = initialNodeX;
              let currentY = initialNodeY;

              const onMouseMove = (moveEvent: MouseEvent) => {
                currentX = initialNodeX + (moveEvent.clientX - startX);
                currentY = initialNodeY + (moveEvent.clientY - startY);
                onUpdateNode({ ...node, x: currentX, y: currentY }, true);
              };
              
              const onMouseUp = () => {
                window.removeEventListener('mousemove', onMouseMove);
                window.removeEventListener('mouseup', onMouseUp);
                if (Math.abs(currentX - initialNodeX) > 1 || Math.abs(currentY - initialNodeY) > 1) {
                  onUpdateNode({ ...node, x: currentX, y: currentY }, false);
                }
              };
              window.addEventListener('mousemove', onMouseMove);
              window.addEventListener('mouseup', onMouseUp);
            }}
            style={{
              left: node.x,
              top: node.y
            }}
            className={`absolute -translate-x-1/2 -translate-y-1/2 p-3 bg-white border-2 rounded-lg shadow-md cursor-grab active:cursor-grabbing transition-all hover:scale-105 z-10
              ${selectedId === node.id ? 'border-blue-500 ring-2 ring-blue-100' : 'border-slate-300'}
              ${isLinking === node.id ? 'border-yellow-500 animate-pulse' : ''}
            `}
          >
            <div className="flex flex-col items-center gap-1 min-w-[50px] pointer-events-none">
              {getNodeIcon(node)}
              <span className="text-xs font-bold text-slate-700 whitespace-nowrap">{node.name}</span>
              <span className="text-[10px] text-slate-400 fira-code">{node.voltage.toFixed(3)} V</span>
            </div>
            <button 
              onClick={(e) => { e.stopPropagation(); setIsLinking(node.id); }}
              className="absolute -right-2 top-1/2 -translate-y-1/2 w-4 h-4 bg-slate-100 border border-slate-300 rounded-full flex items-center justify-center hover:bg-blue-500 hover:text-white transition-colors"
              title="Add connection"
            >
              <Plus size={10} />
            </button>
          </div>
        ))}

      </div>
    </div>
  );
};

export default NetworkCanvas;
