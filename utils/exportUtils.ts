
import { NetworkState, NodeType } from '../types';

/**
 * Converts the current network state to a simplified PSS/E RAW (V33) format.
 * This covers Bus, Load, Generator, and Branch records.
 */
export const convertToPSSERaw = (network: NetworkState): string => {
  const lines: string[] = [];
  const baseMVA = 100.0;

  // Header Records
  lines.push(`0, ${baseMVA}, 33, 0, 0, 60.0  / PSS/E RAW Export from PowerGrid Designer`);
  lines.push(`Network Design Export`);
  lines.push(`Generated on ${new Date().toISOString()}`);

  // 1. BUS DATA
  lines.push(`0 / END OF IC DATA, BEGIN BUS DATA`);
  network.nodes.forEach((node, index) => {
    const busNum = index + 1;
    const type = node.type === NodeType.SLACK ? 3 : node.type === NodeType.PV ? 2 : 1;
    // Format: I, 'NAME', BASKV, IDE, GL, BL, AREA, ZONE, VM, VA, OWNER
    lines.push(`${busNum}, '${node.name.substring(0, 12)}', 138.00, ${type}, 1, 1, 1, ${node.voltage.toFixed(4)}, ${node.angle.toFixed(4)}, 1.1, 0.9`);
  });
  lines.push(`0 / END OF BUS DATA, BEGIN LOAD DATA`);

  // 2. LOAD DATA
  network.nodes.filter(n => n.pLoad > 0 || n.qLoad > 0).forEach((node, index) => {
    const busNum = network.nodes.findIndex(n => n.id === node.id) + 1;
    // Format: I, ID, STATUS, AREA, ZONE, PL, QL, IP, IQ, YP, YQ, OWNER
    lines.push(`${busNum}, '1', 1, 1, 1, ${node.pLoad.toFixed(3)}, ${node.qLoad.toFixed(3)}, 0.0, 0.0, 0.0, 0.0, 1`);
  });
  lines.push(`0 / END OF LOAD DATA, BEGIN FIXED SHUNT DATA`);
  lines.push(`0 / END OF FIXED SHUNT DATA, BEGIN GENERATOR DATA`);

  // 3. GENERATOR DATA
  network.nodes.filter(n => n.type !== NodeType.PQ).forEach((node) => {
    const busNum = network.nodes.findIndex(n => n.id === node.id) + 1;
    // Format: I, ID, PG, QG, QT, QB, VS, IREG, MBASE, ZR, ZX, RT, RX, GTAP, STAT, RMPCT, PT, PB
    lines.push(`${busNum}, '1', ${node.pGen.toFixed(3)}, ${node.qGen.toFixed(3)}, 999.0, -999.0, 1.0, 0, 100.0, 0.0, 0.0, 0.0, 1.0, 1, 100.0, 999.0, -999.0`);
  });
  lines.push(`0 / END OF GENERATOR DATA, BEGIN BRANCH DATA`);

  // 4. BRANCH DATA
  network.edges.forEach((edge) => {
    const fromIdx = network.nodes.findIndex(n => n.id === edge.from) + 1;
    const toIdx = network.nodes.findIndex(n => n.id === edge.to) + 1;
    // Format: I, J, CKTR, R, X, B, RATEA, RATEB, RATEC, GI, BI, GJ, BJ, ST, MET, LEN
    lines.push(`${fromIdx}, ${toIdx}, '1', ${edge.resistance.toFixed(5)}, ${edge.reactance.toFixed(5)}, ${edge.susceptance.toFixed(5)}, ${edge.limit.toFixed(2)}, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1, 1, 0.0`);
  });
  lines.push(`0 / END OF BRANCH DATA, BEGIN TRANSFORMER DATA`);
  lines.push(`0 / END OF TRANSFORMER DATA`);
  lines.push(`Q`);

  return lines.join('\n');
};

/**
 * Trigger a browser download of the provided string data.
 */
export const downloadFile = (data: string, fileName: string, mimeType: string) => {
  const blob = new Blob([data], { type: mimeType });
  const url = window.URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = fileName;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  window.URL.revokeObjectURL(url);
};
