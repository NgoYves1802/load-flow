
import React, { useState, useCallback, useRef, useEffect, useLayoutEffect } from 'react';
import { BusNode, BranchLine, NetworkState, CalculationResult, NodeType, GAParameters, GAResult, SelectionMethod } from './types';
import NetworkCanvas from './components/NetworkCanvas';
import ControlPanel from './components/ControlPanel';
import LibraryPanel, { LibraryItem } from './components/LibraryPanel';
import { calculateLoadFlow } from './services/geminiService';
import { calculateLoadFlowGA } from './services/gaService';
import { IEEE_3_BUS, IEEE_5_BUS, IEEE_14_BUS, RADIAL_FEEDER } from './services/templateService';
import { convertToPSSERaw, downloadFile } from './utils/exportUtils';
import {
  Layout, Database, Terminal, Github, Download, FileJson,
  FileCode, ChevronDown, Undo2, Redo2, BookOpen, FileUp,
  Sparkles, Network, Import, Layers, BarChart3, Activity, Plus, Minus, RotateCcw
} from 'lucide-react';

const App: React.FC = () => {
  const [network, setNetwork] = useState<NetworkState>({
    nodes: [],
    edges: []
  });
  
  const [history, setHistory] = useState<NetworkState[]>([{ nodes: [], edges: [] }]);
  const [historyIndex, setHistoryIndex] = useState(0);

  const [selectedId, setSelectedId] = useState<string | null>(null);
  const [selectedType, setSelectedType] = useState<'node' | 'edge' | null>(null);
  const [isCalculating, setIsCalculating] = useState(false);
  const [results, setResults] = useState<CalculationResult | null>(null);
  const [showExportMenu, setShowExportMenu] = useState(false);
  const [showImportMenu, setShowImportMenu] = useState(false);
  const [isLibraryOpen, setIsLibraryOpen] = useState(true);
  const [gaParams, setGaParams] = useState<GAParameters>({
    populationSize: 100,
    maxGenerations: 300,
    crossoverRate: 0.8,
    mutationRate: 0.05,
    convergenceThreshold: 0.001,
    minVoltage: 0.95,
    maxVoltage: 1.05,
    selectionMethod: SelectionMethod.TOURNAMENT,
    tournamentSize: 5,
    alpha: 1e6, // Increased for delta P violations
    beta: 1e6, // Increased for delta Q violations
    gamma: 1e5 // Increased for Qg violations
  });
  const [activeMainTab, setActiveMainTab] = useState<'schema' | 'results'>('schema');
  const [canvasZoom, setCanvasZoom] = useState(1);

  // Render fitness chart when results change
  useLayoutEffect(() => {
    try {
      if (results?.gaResult && activeMainTab === 'results') {
        const canvas = document.getElementById('fitnessChart') as HTMLCanvasElement;
        if (!canvas) {
          console.warn('Fitness chart canvas not found');
          return;
        }

        const ctx = canvas.getContext('2d');
        if (!ctx) {
          console.warn('Could not get 2D context from canvas');
          return;
        }

        const history = results.gaResult.fitnessHistory;
        if (!Array.isArray(history) || history.length === 0) {
          console.warn('Fitness history is empty or not an array');
          return;
        }

        // Ensure canvas has dimensions
        if (canvas.offsetWidth === 0 || canvas.offsetHeight === 0) {
          console.warn('Canvas has no dimensions');
          return;
        }

        canvas.width = canvas.offsetWidth;
        canvas.height = canvas.offsetHeight;

        const w = canvas.width;
        const h = canvas.height;
        const padding = 50;

        // Validate dimensions
        if (w <= 2 * padding || h <= 2 * padding) {
          console.warn('Canvas too small for chart');
          return;
        }

        ctx.clearRect(0, 0, w, h);

        const maxFit = Math.max(...history);
        const minFit = Math.min(...history);

        // Validate fitness values
        if (!isFinite(maxFit) || !isFinite(minFit)) {
          console.warn('Invalid fitness values');
          return;
        }

        // Draw axes
        ctx.strokeStyle = '#333';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(padding, padding);
        ctx.lineTo(padding, h - padding);
        ctx.lineTo(w - padding, h - padding);
        ctx.stroke();

        // Draw horizontal grid lines
        ctx.strokeStyle = '#e0e0e0';
        ctx.lineWidth = 1;
        for (let i = 1; i < 5; i++) {
          const y = padding + (i / 5) * (h - 2 * padding);
          ctx.beginPath();
          ctx.moveTo(padding, y);
          ctx.lineTo(w - padding, y);
          ctx.stroke();
        }

        // Draw vertical grid lines for generation intervals
        const numVerticalLines = Math.min(10, history.length);
        const verticalInterval = Math.max(1, Math.floor(history.length / numVerticalLines));

        for (let i = 0; i < history.length; i += verticalInterval) {
          const x = padding + (i / (history.length - 1)) * (w - 2 * padding);
          ctx.beginPath();
          ctx.moveTo(x, padding);
          ctx.lineTo(x, h - padding);
          ctx.stroke();
        }

        // Draw curve
        ctx.strokeStyle = '#667eea';
        ctx.lineWidth = 3;
        ctx.beginPath();

        history.forEach((fit, i) => {
          if (typeof fit !== 'number' || !isFinite(fit)) return;

          const x = padding + (i / (history.length - 1)) * (w - 2 * padding);
          const y = h - padding - ((fit - minFit) / (maxFit - minFit || 1)) * (h - 2 * padding);

          if (i === 0) ctx.moveTo(x, y);
          else ctx.lineTo(x, y);
        });

        ctx.stroke();

        // Labels
        ctx.fillStyle = '#333';
        ctx.font = '14px Arial';
        ctx.fillText('Generation', w / 2 - 40, h - 15);
        ctx.save();
        ctx.translate(20, h / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.fillText('Fitness', 0, 0);
        ctx.restore();

        // Generation axis labels (more detailed)
        ctx.font = '11px Arial';
        ctx.fillStyle = '#666';
        ctx.textAlign = 'center';

        // Show generation numbers at regular intervals
        const numLabels = Math.min(10, history.length); // Show up to 10 labels
        const labelInterval = Math.max(1, Math.floor(history.length / numLabels));

        for (let i = 0; i < history.length; i += labelInterval) {
          const x = padding + (i / (history.length - 1)) * (w - 2 * padding);
          const generationNum = i + 1; // Generations are 1-indexed
          ctx.fillText(generationNum.toString(), x, h - padding + 18);
        }

        // Always show the last generation number
        if (history.length > 0) {
          const lastX = padding + (w - 2 * padding);
          ctx.fillText(history.length.toString(), lastX, h - padding + 18);
        }

        ctx.textAlign = 'left'; // Reset text alignment

        // Value labels
        ctx.font = '12px Arial';
        ctx.fillStyle = '#666';
        ctx.fillText(maxFit.toFixed(4), padding - 45, padding + 5);
        ctx.fillText(minFit.toFixed(4), padding - 45, h - padding + 5);
      }
    } catch (error) {
      console.error('Error rendering fitness chart:', error);
    }
  }, [results, activeMainTab]);
  
  const exportMenuRef = useRef<HTMLDivElement>(null);
  const importMenuRef = useRef<HTMLDivElement>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === 'z') {
        if (e.shiftKey) handleRedo(); else handleUndo();
      } else if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === 'y') {
        handleRedo();
      }
    };
    window.addEventListener('keydown', handleKeyDown);
    return () => window.removeEventListener('keydown', handleKeyDown);
  }, [history, historyIndex]);

  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (exportMenuRef.current && !exportMenuRef.current.contains(event.target as Node)) {
        setShowExportMenu(false);
      }
      if (importMenuRef.current && !importMenuRef.current.contains(event.target as Node)) {
        setShowImportMenu(false);
      }
    };
    document.addEventListener('mousedown', handleClickOutside);
    return () => document.removeEventListener('mousedown', handleClickOutside);
  }, []);

  const pushToHistory = (nextState: NetworkState) => {
    const newHistory = history.slice(0, historyIndex + 1);
    newHistory.push(JSON.parse(JSON.stringify(nextState)));
    setHistory(newHistory);
    setHistoryIndex(newHistory.length - 1);
  };

  const handleUndo = () => {
    if (historyIndex > 0) {
      const prevState = history[historyIndex - 1];
      setNetwork(JSON.parse(JSON.stringify(prevState)));
      setHistoryIndex(historyIndex - 1);
    }
  };

  const handleRedo = () => {
    if (historyIndex < history.length - 1) {
      const nextState = history[historyIndex + 1];
      setNetwork(JSON.parse(JSON.stringify(nextState)));
      setHistoryIndex(historyIndex + 1);
    }
  };

  const handleLoadNetwork = (newState: NetworkState) => {
    if (network.nodes.length > 0) {
      if (!confirm("This will replace your current design. Continue?")) return;
    }
    setNetwork(newState);
    pushToHistory(newState);
    setSelectedId(null);
    setResults(null);
    setShowImportMenu(false);
  };

  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (e) => {
      try {
        const content = e.target?.result as string;
        const parsed = JSON.parse(content) as NetworkState;
        if (parsed.nodes && Array.isArray(parsed.nodes)) {
          handleLoadNetwork(parsed);
        } else {
          throw new Error("Invalid format");
        }
      } catch (err) {
        alert("Failed to parse network file. Please ensure it is a valid Design JSON.");
      }
    };
    reader.readAsText(file);
    event.target.value = ''; // Reset
  };

  const handleLibraryItemClick = (item: LibraryItem) => {
    if (item.category === 'Transmission' && (item.id === 'branch-line' || item.id === 'branch-transformer')) {
      alert("To add a Line or Transformer: Click the '+' icon on an existing Bus Bar and then click on the destination Bus Bar.");
      return;
    }

    const spawnX = 300 + (Math.random() * 60 - 30);
    const spawnY = 300 + (Math.random() * 60 - 30);

    const newNode: BusNode = {
      id: `B-${Date.now()}`,
      name: item.defaultData.name || `Bus ${network.nodes.length + 1}`,
      type: item.type,
      category: item.category,
      voltage: item.defaultData.voltage || 1.0,
      angle: 0,
      pGen: item.defaultData.pGen || 0,
      qGen: item.defaultData.qGen || 0,
      pLoad: item.defaultData.pLoad || 0,
      qLoad: item.defaultData.qLoad || 0,
      x: spawnX,
      y: spawnY,
      iconType: item.id
    };

    const next = { ...network, nodes: [...network.nodes, newNode] };
    setNetwork(next);
    pushToHistory(next);
    setSelectedId(newNode.id);
    setSelectedType('node');
  };

  const handleAddNode = (node: BusNode) => {
    const next = { ...network, nodes: [...network.nodes, node] };
    setNetwork(next);
    pushToHistory(next);
  };

  const handleAddEdge = (edge: BranchLine) => {
    const next = { ...network, edges: [...network.edges, edge] };
    setNetwork(next);
    pushToHistory(next);
  };

  const handleUpdateNode = (node: BusNode, isTransient: boolean = false) => {
    setNetwork(prev => {
      const next = {
        ...prev,
        nodes: prev.nodes.map(n => n.id === node.id ? node : n)
      };
      if (!isTransient) pushToHistory(next);
      return next;
    });
  };

  const handleUpdateEdge = (edge: BranchLine) => {
    const next = {
      ...network,
      edges: network.edges.map(e => e.id === edge.id ? edge : e)
    };
    setNetwork(next);
    pushToHistory(next);
  };

  const handleDeleteNode = (id: string) => {
    const next = {
      nodes: network.nodes.filter(n => n.id !== id),
      edges: network.edges.filter(e => e.from !== id && e.to !== id)
    };
    setNetwork(next);
    pushToHistory(next);
    setSelectedId(null);
  };

  const handleDeleteEdge = (id: string) => {
    const next = {
      ...network,
      edges: network.edges.filter(e => e.id !== id)
    };
    setNetwork(next);
    pushToHistory(next);
    setSelectedId(null);
  };

  const runLoadFlow = async () => {
    setIsCalculating(true);
    try {
      // Run GA calculation
      const gaResult = await calculateLoadFlowGA(network, gaParams, (gen, fitness) => {
        // Progress callback - could be used for progress updates
        console.log(`Generation ${gen}, Fitness: ${fitness}`);
      });

      // For MATLAB-inspired approach, generations are not optimized
      // Slack P_gen is determined by power balance, PV Q_gen by voltage control
      const nodes = network.nodes.map((node, index) => ({
        id: node.id,
        voltage: gaResult.solution.voltages[index],
        angle: gaResult.solution.angles[index] || 0,
        pGen: node.pGen - node.pLoad, // Keep original values
        qGen: node.qGen - node.qLoad // Keep original values
      }));

      // Calculate power flows for edges (simplified calculation)
      const edges = network.edges.map(edge => {
        const fromIdx = network.nodes.findIndex(node => node.id === edge.from);
        const toIdx = network.nodes.findIndex(node => node.id === edge.to);

        if (fromIdx === -1 || toIdx === -1) {
          return { id: edge.id, pFlow: 0, qFlow: 0, losses: 0 };
        }

        const V1 = gaResult.solution.voltages[fromIdx];
        const V2 = gaResult.solution.voltages[toIdx];
        const theta1 = gaResult.solution.angles[fromIdx] * Math.PI / 180;
        const theta2 = gaResult.solution.angles[toIdx] * Math.PI / 180;

        // Simplified power flow calculation
        const deltaTheta = theta1 - theta2;
        const denom = edge.resistance * edge.resistance + edge.reactance * edge.reactance;

        // Current magnitude approximation
        const I_real = (V1 * Math.cos(theta1) - V2 * Math.cos(theta2)) / denom;
        const I_imag = (V1 * Math.sin(theta1) - V2 * Math.sin(theta2)) / denom;

        // Power flow from sending end
        const pFlow = V1 * V1 * edge.resistance / denom + V1 * V2 * (edge.resistance * Math.cos(deltaTheta) + edge.reactance * Math.sin(deltaTheta)) / denom;
        const qFlow = -V1 * V1 * (edge.reactance + edge.susceptance) / denom + V1 * V2 * (edge.reactance * Math.cos(deltaTheta) - edge.resistance * Math.sin(deltaTheta)) / denom;

        // Calculate losses
        const losses = (pFlow * pFlow + qFlow * qFlow) * edge.resistance / (V1 * V1);

        return {
          id: edge.id,
          pFlow: Math.abs(pFlow) * 100, // Convert to MW
          qFlow: Math.abs(qFlow) * 100, // Convert to MVAr
          losses: Math.abs(losses) * 100 // Convert to appropriate units
        };
      });

      const solverResults: CalculationResult = {
        nodes,
        edges,
        summary: `GA converged ${gaResult.converged ? 'successfully' : 'with limited success'}. Final fitness: ${gaResult.fitness.toFixed(6)} after ${gaResult.generations} generations.`,
        gaResult
      };

      setResults(solverResults);

      // Update network with results
      const updatedNodes = network.nodes.map(n => {
        const res = solverResults.nodes.find(r => r.id === n.id);
        return res ? { ...n, voltage: res.voltage, angle: res.angle } : n;
      });

      const next = { nodes: updatedNodes, edges: network.edges };
      setNetwork(next);
      pushToHistory(next);
    } catch (error) {
      console.error(error);
      alert("Calculation failed. Please check network topology and GA parameters.");
    } finally {
      setIsCalculating(false);
    }
  };

  const exportAsJSON = () => {
    const data = JSON.stringify(network, null, 2);
    downloadFile(data, 'network-design.json', 'application/json');
    setShowExportMenu(false);
  };

  const exportAsRAW = () => {
    const data = convertToPSSERaw(network);
    downloadFile(data, 'network-design.raw', 'text/plain');
    setShowExportMenu(false);
  };

  return (
    <div className="flex flex-col h-screen bg-slate-100 overflow-hidden">
      <input 
        type="file" 
        ref={fileInputRef} 
        onChange={handleFileUpload} 
        accept=".json" 
        className="hidden" 
      />
      
      <header className="h-14 bg-white border-b border-slate-200 flex items-center justify-between px-6 shadow-sm z-20">
        <div className="flex items-center gap-3">
          <div className="w-8 h-8 bg-blue-600 rounded-lg flex items-center justify-center text-white shadow-md">
            <Layout size={18} />
          </div>
          <div>
            <h1 className="text-sm font-bold text-slate-800 leading-tight">PowerGrid Designer</h1>
            <p className="text-[10px] text-slate-400 font-medium uppercase tracking-tight">System Layout & Load Flow Solver</p>
          </div>
        </div>

        <div className="flex items-center gap-6">
          <div className="flex items-center gap-2">
            <button 
              onClick={handleUndo}
              disabled={historyIndex === 0}
              className="p-1.5 text-slate-500 hover:text-blue-600 hover:bg-slate-50 rounded-md disabled:opacity-30 transition-all"
              title="Undo"
            >
              <Undo2 size={18} />
            </button>
            <button 
              onClick={handleRedo}
              disabled={historyIndex >= history.length - 1}
              className="p-1.5 text-slate-500 hover:text-blue-600 hover:bg-slate-50 rounded-md disabled:opacity-30 transition-all"
              title="Redo"
            >
              <Redo2 size={18} />
            </button>
          </div>

          <div className="w-px h-6 bg-slate-200"></div>

          <nav className="flex items-center gap-4">
            <button 
              onClick={() => setIsLibraryOpen(!isLibraryOpen)}
              className={`flex items-center gap-1.5 text-xs font-semibold transition-colors px-3 py-1.5 rounded-md border ${isLibraryOpen ? 'bg-blue-50 text-blue-600 border-blue-200' : 'bg-slate-50 text-slate-700 border-slate-200 hover:bg-slate-100'}`}
            >
              <BookOpen size={14} /> Library
            </button>

            {/* Import Menu */}
            <div className="relative" ref={importMenuRef}>
              <button 
                onClick={() => setShowImportMenu(!showImportMenu)}
                className="text-slate-700 hover:text-blue-600 flex items-center gap-1.5 text-xs font-semibold transition-colors px-3 py-1.5 bg-slate-50 hover:bg-blue-50 rounded-md border border-slate-200"
              >
                <Import size={14} /> Import <ChevronDown size={12} className={`transition-transform duration-200 ${showImportMenu ? 'rotate-180' : ''}`} />
              </button>
              
              {showImportMenu && (
                <div className="absolute top-full right-0 mt-1 w-56 bg-white border border-slate-200 rounded-lg shadow-xl py-1 z-30 animate-in fade-in slide-in-from-top-2">
                  <button onClick={() => fileInputRef.current?.click()} className="w-full px-4 py-2 text-left text-xs text-slate-600 hover:bg-slate-50 flex items-center gap-3 transition-colors">
                    <FileUp size={14} className="text-blue-500" />
                    <span>Load Design File (.json)</span>
                  </button>
                  <div className="border-t border-slate-100 my-1"></div>
                  <div className="px-4 py-1 text-[10px] font-bold text-slate-400 uppercase tracking-wider">Example Templates</div>
                  <button onClick={() => handleLoadNetwork(IEEE_3_BUS)} className="w-full px-4 py-2 text-left text-xs text-slate-600 hover:bg-slate-50 flex items-center gap-3 transition-colors">
                    <Sparkles size={14} className="text-amber-500" />
                    <span>IEEE 3-Bus System</span>
                  </button>
                  <button onClick={() => handleLoadNetwork(IEEE_5_BUS)} className="w-full px-4 py-2 text-left text-xs text-slate-600 hover:bg-slate-50 flex items-center gap-3 transition-colors">
                    <Network size={14} className="text-emerald-500" />
                    <span>IEEE 5-Bus System</span>
                  </button>
                  <button onClick={() => handleLoadNetwork(IEEE_14_BUS)} className="w-full px-4 py-2 text-left text-xs text-slate-600 hover:bg-slate-50 flex items-center gap-3 transition-colors">
                    <Activity size={14} className="text-purple-500" />
                    <span>IEEE 14-Bus System</span>
                  </button>
                  <button onClick={() => handleLoadNetwork(RADIAL_FEEDER)} className="w-full px-4 py-2 text-left text-xs text-slate-600 hover:bg-slate-50 flex items-center gap-3 transition-colors">
                    <Database size={14} className="text-slate-500" />
                    <span>Radial Feeder Example</span>
                  </button>
                </div>
              )}
            </div>

            {/* Export Menu */}
            <div className="relative" ref={exportMenuRef}>
              <button 
                onClick={() => setShowExportMenu(!showExportMenu)}
                className="text-slate-700 hover:text-blue-600 flex items-center gap-1.5 text-xs font-semibold transition-colors px-3 py-1.5 bg-slate-50 hover:bg-blue-50 rounded-md border border-slate-200"
              >
                <Download size={14} /> Export <ChevronDown size={12} className={`transition-transform duration-200 ${showExportMenu ? 'rotate-180' : ''}`} />
              </button>
              
              {showExportMenu && (
                <div className="absolute top-full right-0 mt-1 w-48 bg-white border border-slate-200 rounded-lg shadow-xl py-1 z-30 animate-in fade-in slide-in-from-top-2">
                  <button onClick={exportAsJSON} className="w-full px-4 py-2 text-left text-xs text-slate-600 hover:bg-slate-50 flex items-center gap-3 transition-colors">
                    <FileJson size={14} className="text-amber-500" />
                    <span>JSON Config (.json)</span>
                  </button>
                  <button onClick={exportAsRAW} className="w-full px-4 py-2 text-left text-xs text-slate-600 hover:bg-slate-50 flex items-center gap-3 transition-colors">
                    <FileCode size={14} className="text-blue-500" />
                    <span>PSS/E RAW V33 (.raw)</span>
                  </button>
                </div>
              )}
            </div>
            
            <button className="text-slate-500 hover:text-blue-600 flex items-center gap-1.5 text-xs font-semibold transition-colors">
              <Database size={14} /> Project
            </button>
            <button className="text-slate-500 hover:text-blue-600 flex items-center gap-1.5 text-xs font-semibold transition-colors">
              <Terminal size={14} /> Console
            </button>
          </nav>
          <div className="w-px h-6 bg-slate-200"></div>
          <button className="p-2 text-slate-400 hover:text-slate-600 transition-colors">
            <Github size={20} />
          </button>
        </div>
      </header>

      <main className="flex-1 flex overflow-hidden p-4 gap-4">
        {isLibraryOpen && <LibraryPanel onItemClick={handleLibraryItemClick} />}

        <div className="flex-1 flex flex-col bg-white border border-slate-200 rounded-lg overflow-hidden">
          {/* Main Drawing Area Tabs */}
          <div className="flex border-b border-slate-200 bg-slate-50">
            <button
              onClick={() => setActiveMainTab('schema')}
              className={`flex-1 py-3 px-4 text-sm font-medium transition-colors flex items-center justify-center gap-2 ${
                activeMainTab === 'schema'
                  ? 'text-blue-600 border-b-2 border-blue-600 bg-white'
                  : 'text-slate-500 hover:text-slate-700 hover:bg-slate-100'
              }`}
            >
              <Layers size={16} />
              Schema
            </button>
            <button
              onClick={() => setActiveMainTab('results')}
              className={`flex-1 py-3 px-4 text-sm font-medium transition-colors flex items-center justify-center gap-2 ${
                activeMainTab === 'results'
                  ? 'text-blue-600 border-b-2 border-blue-600 bg-white'
                  : 'text-slate-500 hover:text-slate-700 hover:bg-slate-100'
              }`}
            >
              <BarChart3 size={16} />
              Results
            </button>
          </div>

          {/* Tab Content */}
          <div className="flex-1 overflow-hidden flex flex-col">
            {activeMainTab === 'schema' && (
              <>
                {/* Drawing Area Header */}
                <div className="px-4 py-2 bg-slate-50 border-b border-slate-200 flex items-center justify-between">
                  <div className="flex items-center gap-2">
                    <Layers size={16} className="text-blue-600" />
                    <span className="text-sm font-medium text-slate-700">Network Design Area</span>
                    <span className="text-xs text-slate-500">• Click components from the library to add them here</span>
                  </div>
                  <div className="flex items-center gap-3">
                    <span className="text-xs text-slate-500">Zoom:</span>
                    <div className="flex items-center gap-1">
                      <button
                        onClick={() => setCanvasZoom(prev => Math.min(prev * 1.2, 5))}
                        className="w-8 h-8 bg-white border border-slate-300 rounded flex items-center justify-center hover:bg-blue-50 hover:border-blue-400 transition-all"
                        title="Zoom In"
                      >
                        <Plus size={14} className="text-slate-600" />
                      </button>
                      <button
                        onClick={() => setCanvasZoom(prev => Math.max(prev / 1.2, 0.1))}
                        className="w-8 h-8 bg-white border border-slate-300 rounded flex items-center justify-center hover:bg-blue-50 hover:border-blue-400 transition-all"
                        title="Zoom Out"
                      >
                        <Minus size={14} className="text-slate-600" />
                      </button>
                      <button
                        onClick={() => setCanvasZoom(1)}
                        className="w-8 h-8 bg-white border border-slate-300 rounded flex items-center justify-center hover:bg-blue-50 hover:border-blue-400 transition-all"
                        title="Reset Zoom"
                      >
                        <RotateCcw size={14} className="text-slate-600" />
                      </button>
                      <div className="w-12 h-6 bg-white border border-slate-300 rounded flex items-center justify-center text-xs font-medium text-slate-600 ml-2">
                        {Math.round(canvasZoom * 100)}%
                      </div>
                    </div>
                  </div>
                </div>

                {/* Drawing Canvas */}
                <div className="flex-1 relative">
                  <NetworkCanvas
                    nodes={network.nodes}
                    edges={network.edges}
                    onAddNode={handleAddNode}
                    onAddEdge={handleAddEdge}
                    onUpdateNode={handleUpdateNode}
                    onUpdateEdge={handleUpdateEdge}
                    onDeleteNode={handleDeleteNode}
                    onDeleteEdge={handleDeleteEdge}
                    selectedId={selectedId}
                    onSelect={(id, type) => { setSelectedId(id); setSelectedType(type); }}
                    zoom={canvasZoom}
                    onZoomChange={setCanvasZoom}
                  />

                  {/* Instructions overlay for empty canvas */}
                  {network.nodes.length === 0 && (
                    <div className="absolute top-4 right-4 pointer-events-none">
                      <div className="bg-blue-50 border border-blue-200 rounded-lg p-4 max-w-xs shadow-sm">
                        <div className="flex items-center gap-2 mb-2">
                          <Layers size={16} className="text-blue-600" />
                          <span className="text-sm font-medium text-blue-700">Drawing Area</span>
                        </div>
                        <p className="text-xs text-blue-600">
                          Click on components from the Library panel to add them here. Drag to connect components.
                        </p>
                      </div>
                    </div>
                  )}
                </div>
              </>
            )}

            {activeMainTab === 'results' && (
              <div className="h-full overflow-y-auto p-6 space-y-6">
                {!results && (
                  <div className="flex flex-col items-center justify-center h-64 text-slate-400 text-center">
                    <BarChart3 size={40} className="mb-4 opacity-20" />
                    <p className="text-sm">Run the solver to see GA convergence results and fitness history.</p>
                  </div>
                )}

                {results && (
                  <div className="space-y-6">
                    {/* GA Performance Summary */}
                    {results.gaResult && (
                      <div className="p-6 bg-blue-50 border border-blue-200 rounded-lg">
                        <h3 className="text-lg font-bold text-blue-800 mb-4">GA Performance Summary</h3>
                        <div className="grid grid-cols-2 gap-6 text-sm">
                          <div>
                            <span className="text-slate-600">Final Fitness:</span>
                            <span className="font-mono text-blue-700 ml-2 text-lg">
                              {typeof results.gaResult.fitness === 'number' && isFinite(results.gaResult.fitness)
                                ? results.gaResult.fitness.toFixed(6)
                                : 'N/A'}
                            </span>
                          </div>
                          <div>
                            <span className="text-slate-600">Generations:</span>
                            <span className="font-mono text-blue-700 ml-2 text-lg">
                              {typeof results.gaResult.generations === 'number' && isFinite(results.gaResult.generations)
                                ? results.gaResult.generations
                                : 'N/A'}
                            </span>
                          </div>
                          <div>
                            <span className="text-slate-600">Converged:</span>
                            <span className={`font-bold ml-2 text-lg ${results.gaResult.converged ? 'text-green-600' : 'text-orange-600'}`}>
                              {typeof results.gaResult.converged === 'boolean'
                                ? (results.gaResult.converged ? 'Yes' : 'No')
                                : 'Unknown'}
                            </span>
                          </div>
                        </div>
                      </div>
                    )}

                    {/* Analysis Report */}
                    <div className="p-4 bg-blue-50 border border-blue-100 rounded-lg">
                      <h3 className="text-lg font-bold text-blue-600 flex items-center gap-2 mb-3">
                        <Activity size={20} /> Analysis Report
                      </h3>
                      <div className="text-sm leading-relaxed text-slate-700 italic">
                        "{results.summary}"
                      </div>
                    </div>

                    {/* Calculation Results */}
                    <div className="space-y-6">
                      <h3 className="text-xl font-bold text-slate-800">Calculation Results</h3>

                      {/* Bus Results */}
                      <div className="bg-white border border-slate-200 rounded-lg overflow-hidden shadow-sm">
                        <div className="px-6 py-4 bg-slate-50 border-b border-slate-200">
                          <h4 className="text-lg font-bold text-slate-700">Bus Results</h4>
                        </div>
                        <div className="overflow-x-auto">
                          <table className="w-full">
                            <thead className="bg-slate-100">
                              <tr>
                                
                                <th className="px-6 py-3 text-left text-sm font-bold text-slate-600">Bus ID</th>
                                <th className="px-6 py-3 text-left text-sm font-bold text-slate-600">Type</th>
                                <th className="px-6 py-3 text-left text-sm font-bold text-slate-600">V (p.u.)</th>
                                <th className="px-6 py-3 text-left text-sm font-bold text-slate-600">Angle (°)</th>
                                <th className="px-6 py-3 text-left text-sm font-bold text-slate-600">P Gen (MW)</th>
                                <th className="px-6 py-3 text-left text-sm font-bold text-slate-600">Q Gen (MVAr)</th>
                                <th className="px-6 py-3 text-left text-sm font-bold text-slate-600">P Cal (MW)</th>
                                <th className="px-6 py-3 text-left text-sm font-bold text-slate-600">Q Cal (MVAr)</th>
                              </tr>
                            </thead>
                        <tbody>
                          {Array.isArray(results.nodes) && results.nodes.map((nodeResult) => {
                            const node = network.nodes.find(n => n.id === nodeResult.id);
                            return (
                              <tr key={nodeResult.id || Math.random()} className="border-t border-slate-100 hover:bg-slate-50">
                                <td className="px-6 py-3 text-sm font-mono text-slate-700">{nodeResult.id || 'Unknown'}</td>
                                <td className="px-6 py-3 text-sm text-slate-700">{node?.type || 'Unknown'}</td>
                                <td className="px-6 py-3 text-sm font-mono text-blue-600 font-bold">
                                  {typeof nodeResult.voltage === 'number' && isFinite(nodeResult.voltage)
                                    ? nodeResult.voltage.toFixed(4)
                                    : 'N/A'}
                                </td>
                                <td className="px-6 py-3 text-sm font-mono text-green-600 font-bold">
                                  {typeof nodeResult.angle === 'number' && isFinite(nodeResult.angle)
                                    ? nodeResult.angle.toFixed(2)
                                    : 'N/A'}
                                </td>
                                <td className="px-6 py-3 text-sm font-mono text-slate-700">
                                  {typeof nodeResult.pGen === 'number' && isFinite(nodeResult.pGen)
                                    ? nodeResult.pGen
                                    : 'N/A'}
                                </td>
                                <td className="px-6 py-3 text-sm font-mono text-slate-700">
                                  {typeof nodeResult.qGen === 'number' && isFinite(nodeResult.qGen)
                                    ? nodeResult.qGen
                                    : 'N/A'}
                                </td>
                                <td className="px-6 py-3 text-sm font-mono text-purple-600 font-bold">
                                  {(() => {
                                    // Calculate active power injection for this bus
                                    const busIdx = network.nodes.findIndex(n => n.id === nodeResult.id);
                                    if (busIdx === -1 || !results.gaResult) return 'N/A';

                                    try {
                                      const V = results.gaResult.solution.voltages;
                                      const theta = results.gaResult.solution.angles.map(a => a * Math.PI / 180);

                                      let P = 0;
                                      const n = network.nodes.length;

                                      // Reconstruct Y-bus for this calculation
                                      const Ybus: any[][] = [];
                                      for (let i = 0; i < n; i++) {
                                        Ybus[i] = [];
                                        for (let j = 0; j < n; j++) {
                                          Ybus[i][j] = { real: 0, imag: 0 };
                                        }
                                      }

                                      // Build Y-bus matrix
                                      network.edges.forEach(edge => {
                                        const fromIdx = network.nodes.findIndex(node => node.id === edge.from);
                                        const toIdx = network.nodes.findIndex(node => node.id === edge.to);
                                        if (fromIdx === -1 || toIdx === -1) return;

                                        const Z_mag_sq = edge.resistance * edge.resistance + edge.reactance * edge.reactance;
                                        if (Z_mag_sq === 0) return;

                                        const y_real = edge.resistance / Z_mag_sq;
                                        const y_imag = -edge.reactance / Z_mag_sq;
                                        const y_shunt = edge.susceptance / 2;

                                        Ybus[fromIdx][fromIdx].real += y_real;
                                        Ybus[fromIdx][fromIdx].imag += y_imag - y_shunt;
                                        Ybus[toIdx][toIdx].real += y_real;
                                        Ybus[toIdx][toIdx].imag += y_imag - y_shunt;
                                        Ybus[fromIdx][toIdx].real -= y_real;
                                        Ybus[fromIdx][toIdx].imag -= y_imag;
                                        Ybus[toIdx][fromIdx].real -= y_real;
                                        Ybus[toIdx][fromIdx].imag -= y_imag;
                                      });

                                      // Calculate power injection
                                      for (let j = 0; j < n; j++) {
                                        const thetaIJ = theta[busIdx] - theta[j];
                                        const G = Ybus[busIdx][j].real;
                                        const B = Ybus[busIdx][j].imag;

                                        P += V[j] * (G * Math.cos(thetaIJ) + B * Math.sin(thetaIJ));
                                      }

                                      P *= V[busIdx];
                                      P *= 100; // Convert to MW
                                      return P.toFixed(2);
                                    } catch (error) {
                                      return 'Error';
                                    }
                                  })()}
                                </td>
                                <td className="px-6 py-3 text-sm font-mono text-orange-600 font-bold">
                                  {(() => {
                                    // Calculate reactive power injection for this bus
                                    const busIdx = network.nodes.findIndex(n => n.id === nodeResult.id);
                                    if (busIdx === -1 || !results.gaResult) return 'N/A';

                                    try {
                                      const V = results.gaResult.solution.voltages;
                                      const theta = results.gaResult.solution.angles.map(a => a * Math.PI / 180);

                                      let Q = 0;
                                      const n = network.nodes.length;

                                      // Reconstruct Y-bus for this calculation (same as above)
                                      const Ybus: any[][] = [];
                                      for (let i = 0; i < n; i++) {
                                        Ybus[i] = [];
                                        for (let j = 0; j < n; j++) {
                                          Ybus[i][j] = { real: 0, imag: 0 };
                                        }
                                      }

                                      // Build Y-bus matrix
                                      network.edges.forEach(edge => {
                                        const fromIdx = network.nodes.findIndex(node => node.id === edge.from);
                                        const toIdx = network.nodes.findIndex(node => node.id === edge.to);
                                        if (fromIdx === -1 || toIdx === -1) return;

                                        const Z_mag_sq = edge.resistance * edge.resistance + edge.reactance * edge.reactance;
                                        if (Z_mag_sq === 0) return;

                                        const y_real = edge.resistance / Z_mag_sq;
                                        const y_imag = -edge.reactance / Z_mag_sq;
                                        const y_shunt = edge.susceptance / 2;

                                        Ybus[fromIdx][fromIdx].real += y_real;
                                        Ybus[fromIdx][fromIdx].imag += y_imag - y_shunt;
                                        Ybus[toIdx][toIdx].real += y_real;
                                        Ybus[toIdx][toIdx].imag += y_imag - y_shunt;
                                        Ybus[fromIdx][toIdx].real -= y_real;
                                        Ybus[fromIdx][toIdx].imag -= y_imag;
                                        Ybus[toIdx][fromIdx].real -= y_real;
                                        Ybus[toIdx][fromIdx].imag -= y_imag;
                                      });

                                      // Calculate reactive power injection
                                      for (let j = 0; j < n; j++) {
                                        const thetaIJ = theta[busIdx] - theta[j];
                                        const G = Ybus[busIdx][j].real;
                                        const B = Ybus[busIdx][j].imag;

                                        Q += V[j] * (G * Math.sin(thetaIJ) - B * Math.cos(thetaIJ));
                                      }

                                      Q *= V[busIdx];
                                      Q *= 100; // Convert to MVAr
                                      return Q.toFixed(2);
                                    } catch (error) {
                                      return 'Error';
                                    }
                                  })()}
                                </td>
                              </tr>
                            );
                          })}
                        </tbody>
                          </table>
                        </div>
                      </div>

                      {/* Line Results */}

                    </div>

                    {/* Fitness History Chart */}
                    {results.gaResult && Array.isArray(results.gaResult.fitnessHistory) && results.gaResult.fitnessHistory.length > 0 && (
                      <div className="space-y-4">
                        <h3 className="text-xl font-bold text-slate-800">Fitness Convergence</h3>
                        <div className="h-[500px] bg-slate-50 border border-slate-200 rounded-lg p-6 shadow-sm">
                          <canvas id="fitnessChart" className="w-full h-full"></canvas>
                        </div>
                      </div>
                    )}
                  </div>
                )}
              </div>
            )}
          </div>
        </div>

        <ControlPanel
          selectedId={selectedId}
          selectedType={selectedType}
          nodes={network.nodes}
          edges={network.edges}
          onUpdateNode={handleUpdateNode}
          onUpdateEdge={handleUpdateEdge}
          onDeleteNode={handleDeleteNode}
          onDeleteEdge={handleDeleteEdge}
          onRunLoadFlow={runLoadFlow}
          isCalculating={isCalculating}
          results={results}
          gaParams={gaParams}
          onUpdateGAParams={setGaParams}
        />
      </main>

      <footer className="h-8 bg-white border-t border-slate-200 px-4 flex items-center justify-between text-[10px] text-slate-400 font-medium">
        <div className="flex items-center gap-4">
          <span className="flex items-center gap-1.5"><span className="w-1.5 h-1.5 rounded-full bg-green-500"></span> Engine Ready</span>
          <span>Buses: {network.nodes.length}</span>
          <span>Lines: {network.edges.length}</span>
          <div className="w-px h-3 bg-slate-200"></div>
          <span>Rev: {historyIndex + 1}/{history.length}</span>
        </div>
        <div className="flex items-center gap-4">
          <span>Base S: 100 MVA</span>
          <span className="text-blue-500 font-bold uppercase tracking-widest">GEMINI AI SOLVER</span>
        </div>
      </footer>
    </div>
  );
};

export default App;
