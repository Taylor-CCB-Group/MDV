import { IAnyType } from "mobx-state-tree";
import jsonSchemaToMST from "./jsonSchemaToMST";
import React from "react";

/**
 * Represents a registered chart type, including its MST model and React component.
 */
type ChartTypeEntry = {
  /** The MobX-State-Tree model generated from the chart's JSON schema */
  chartModel: IAnyType;
  /** The React component used to render the chart */
  ChartComponent: React.ComponentType<any>;
};

/**
 * Stores all registered chart types by name.
 */
const allChartTypes: Record<string, ChartTypeEntry> = {};

/**
 * Registers a chart type by name, associating it with a React component and a JSON schema.
 * The schema is converted to an MST model for state management.
 *
 * @param name - The unique name of the chart type.
 * @param options - An object containing the React component and the JSON schema for the chart.
 */
function registerChart(
  name: string,
  { ChartComponent, schema }: { ChartComponent: React.ComponentType<any>; schema: any }
) {
  allChartTypes[name] = {
    chartModel: jsonSchemaToMST(schema),
    ChartComponent,
  };
}

/**
 * Retrieves the registered chart type entry for a given name.
 *
 * @param name - The name of the chart type to retrieve.
 * @returns The ChartTypeEntry if found, otherwise undefined.
 */
function getChart(
  name: string
): ChartTypeEntry | undefined {
  return allChartTypes[name];
}

export { registerChart, getChart };