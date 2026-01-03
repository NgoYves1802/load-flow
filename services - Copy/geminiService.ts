
import { GoogleGenAI, Type } from "@google/genai";
import { NetworkState, CalculationResult } from "../types";

const ai = new GoogleGenAI({ apiKey: process.env.API_KEY });

export const calculateLoadFlow = async (network: NetworkState): Promise<CalculationResult> => {
  const model = 'gemini-3-pro-preview';
  
  const prompt = `Perform a standard power system Load Flow calculation for the following electrical network. 
  Nodes: ${JSON.stringify(network.nodes)}
  Edges: ${JSON.stringify(network.edges)}
  
  Assume a base of 100 MVA. 
  Provide the calculated voltages (p.u.), angles (degrees), power flows (MW/MVAr), and losses (MW).
  Also provide a brief engineering summary of the system stability and any violations (overloads or voltage drops).`;

  const response = await ai.models.generateContent({
    model,
    contents: prompt,
    config: {
      responseMimeType: "application/json",
      responseSchema: {
        type: Type.OBJECT,
        properties: {
          nodes: {
            type: Type.ARRAY,
            items: {
              type: Type.OBJECT,
              properties: {
                id: { type: Type.STRING },
                voltage: { type: Type.NUMBER },
                angle: { type: Type.NUMBER },
                pGen: { type: Type.NUMBER },
                qGen: { type: Type.NUMBER }
              },
              required: ["id", "voltage", "angle"]
            }
          },
          edges: {
            type: Type.ARRAY,
            items: {
              type: Type.OBJECT,
              properties: {
                id: { type: Type.STRING },
                pFlow: { type: Type.NUMBER },
                qFlow: { type: Type.NUMBER },
                losses: { type: Type.NUMBER }
              },
              required: ["id", "pFlow", "qFlow", "losses"]
            }
          },
          summary: { type: Type.STRING }
        },
        required: ["nodes", "edges", "summary"]
      }
    }
  });

  try {
    return JSON.parse(response.text);
  } catch (error) {
    console.error("Failed to parse calculation results:", error);
    throw new Error("Invalid response format from solver engine.");
  }
};
