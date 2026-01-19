import type DataStore from "@/datastore/DataStore";
import type { Gate } from "./types";
import { action, computed, makeObservable, observable } from "mobx";
import type { LoadedDataColumn } from "@/charts/charts";
import { extractCoords, isPointInGate } from "./gateUtils";
import { loadColumn } from "@/dataloaders/DataLoaderUtil";

const GATES_COLUMN_NAME = "__gates__";
//todo: Add update and delete operations
/**
 * Class which manages the gating operations
 * Creates a gates column if it doesn't exist, loads the gates column if it's not loaded
 * Adds new gates and assigns it to the respective cells in the column
 */
export class GateManager {
    gates: Map<string, Gate> = new Map();

    private dataStore: DataStore;
    private gateColumn: LoadedDataColumn<"multitext"> | null = null;

    constructor(dataStore: DataStore) {
        this.dataStore = dataStore;

        makeObservable(this, {
            gates: observable,
            addGate: action,
            gatesArray: computed,
        });

        this.initializeGateColumn();
        this.loadGatesFromConfig();
    }

    /**
     * Ensure gate column's data is loaded, load the column data explicity if not loaded
     */
    private async ensureGateColumnDataLoaded() {
        if (!this.gateColumn) return;

        // Check if data is already loaded
        if (this.gateColumn.data) {
            return;
        }

        // Check if column exists in columnIndex
        if (!this.dataStore.columnIndex[GATES_COLUMN_NAME]) {
            return;
        }

        // Explicitly load the column data
        try {
            await loadColumn(this.dataStore.name, GATES_COLUMN_NAME);
            this.gateColumn = this.dataStore.columnIndex[GATES_COLUMN_NAME] as LoadedDataColumn<"multitext">;
        } catch (error) {
            console.error("Failed to load gates column data:", error);
        }
    }

    /**
     * Fetch the gates column if it exists, ensure the column data is loaded
     * Else create a new column and add the column to datastore
     */
    private initializeGateColumn() {
        const existingCol = this.dataStore.columnIndex[GATES_COLUMN_NAME];
        if (existingCol && existingCol.datatype === "multitext") {
            this.gateColumn = existingCol as LoadedDataColumn<"multitext">;

            // Ensure the column has data loaded
            if (!this.gateColumn.data) {
                console.warn("Gates column exists but data is not loaded yet");
            }

            //! Check if this is right
            // Ensure stringLength is set if missing
            if (typeof this.gateColumn.stringLength !== "number") {
                // @ts-expect-error: Forcibly assign stringLength if not set
                this.gateColumn.stringLength = 10 as any;
            }

            // Ensure values array exists
            if (!Array.isArray(this.gateColumn.values)) {
                this.gateColumn.values = [];
            }

            if (!this.gateColumn.data) {
                // Load data asynchronously
                this.ensureGateColumnDataLoaded().catch(console.error);
            }
        } else {
            const column = {
                name: "Gates",
                field: GATES_COLUMN_NAME,
                datatype: "multitext" as const,
                editable: false, // as of now
                delimiter: "," as const,
                values: [""] as string[],
                stringLength: 10,
            };

            const data = new SharedArrayBuffer(this.dataStore.size * column.stringLength * 2);
            const dataArray = new Uint16Array(data);

            dataArray.fill(65535);

            this.dataStore.addColumn(column, data, true);

            this.gateColumn = this.dataStore.columnIndex[GATES_COLUMN_NAME] as LoadedDataColumn<"multitext">;
        }
    }

    /**
     * Get the gates saved in the dataStore.config or dataStore
     * Ensure the gate column's data is loaded
     * Update the gates column
     */
    private async loadGatesFromConfig() {
        const savedGates = this.dataStore.config.gates || this.dataStore.gates;
        if (Array.isArray(savedGates)) {
            for (const gate of savedGates) {
                const clonedGate = {
                    ...gate,
                    geometry: JSON.parse(JSON.stringify(gate.geometry)),
                };
                this.gates.set(gate.id, clonedGate);
            }
        }

        await this.updateGatesColumnWhenReady();
    }

    /**
     * Utility function to get or add the value index from gate column for the given gate name
     */
    private getOrAddValueIndex(gateName: string): number {
        if (!this.gateColumn) {
            throw new Error("Gate column not initialized");
        }

        let index = this.gateColumn.values.indexOf(gateName);
        if (index === -1) {
            this.gateColumn.values.push(gateName);
            if (this.gateColumn.values.length > 65536) {
                throw new Error(`Gates column exceeded 65536 values when adding '${gateName}'`);
            }
            index = this.gateColumn.values.length - 1;
        }
        return index;
    }

    /**
     * Get the gate names for a given cell
     *
     * @param cellIndex - index of the cell in column's data array
     * @returns array of strings with the gate names
     */
    private getGateNamesForCell(cellIndex: number): string[] {
        if (!this.gateColumn || !this.gateColumn.stringLength) return [];

        const baseIndex = cellIndex * this.gateColumn.stringLength;
        const gateNames: string[] = [];

        for (let i = 0; i < this.gateColumn.stringLength; i++) {
            const valueIndex = this.gateColumn.data[baseIndex + i];
            if (valueIndex === 65535) break;
            if (valueIndex >= 0 && valueIndex < this.gateColumn.values.length) {
                gateNames.push(this.gateColumn.values[valueIndex]);
            }
        }

        return gateNames;
    }

    /**
     * Set the gate names for a given cell with the updated gate names
     *
     * @param cellIndex - index of the cell in the data array
     * @param gateNames - updated gates names associated with the cell
     */
    private setGateNamesForCell(cellIndex: number, gateNames: string[]) {
        if (!this.gateColumn || !this.gateColumn.stringLength) return;

        const baseIndex = cellIndex * this.gateColumn.stringLength;
        const maxValues = Math.min(gateNames.length, this.gateColumn.stringLength);

        for (let i = 0; i < this.gateColumn.stringLength; i++) {
            this.gateColumn.data[baseIndex + i] = 65535;
        }

        const sortedNames = [...gateNames].sort();

        for (let i = 0; i < maxValues; i++) {
            const valueIndex = this.getOrAddValueIndex(sortedNames[i]);
            this.gateColumn.data[baseIndex + i] = valueIndex;
        }
    }

    /**
     * Ensure the gate column data is loaded
     * Set the gate names as empty for all cells
     * Recompute gates for the cells
     */
    private async updateGateColumn() {
        if (!this.gateColumn) {
            return;
        }

        // Ensure data is loaded
        if (!this.gateColumn.data) {
            await this.ensureGateColumnDataLoaded();
            // If still not loaded, can't update yet
            if (!this.gateColumn.data) {
                console.warn("Cannot update gates column: data not loaded");
                return;
            }
        }

        for (let i = 0; i < this.dataStore.size; i++) {
            this.setGateNamesForCell(i, []);
        }

        let anyColumnsLoaded = false;
        for (const gate of this.gates.values()) {
            const [xField, yField] = gate.columns;
            const xCol = this.dataStore.columnIndex[xField];
            const yCol = this.dataStore.columnIndex[yField];

            // Only update if columns exist and are loaded
            if (xCol?.data && yCol?.data) {
                anyColumnsLoaded = true;
                this.updateCellsWithGate(gate, true);
            }
        }

        // Only notify something was updated
        if (anyColumnsLoaded) {
            this.dataStore.dataChanged([GATES_COLUMN_NAME]);
        }
    }

    /**
     * Update the cells in the gate column with the given Gate
     * Add the gate to the cells inside the selection if 'add' param is true
     * Remove the gate from the cells if 'add' is false, if the cells had the gate stored
     */
    private updateCellsWithGate(gate: Gate, add: boolean) {
        if (!this.gateColumn) return;

        const [xField, yField] = gate.columns;

        const xCol = this.dataStore.columnIndex[xField];
        const yCol = this.dataStore.columnIndex[yField];

        if (!xCol || !yCol || !xCol?.data || !yCol?.data) return;

        const polygonCoords = extractCoords(gate.geometry);
        if (polygonCoords.length === 0) return;

        for (let i = 0; i < this.dataStore.size; i++) {
            const x = xCol.data[i];
            const y = yCol.data[i];

            if (!Number.isFinite(x) || !Number.isFinite(y)) continue;

            const gateNames = this.getGateNamesForCell(i);

            let newGateNames: string[];

            if (add) {
                const isInside = isPointInGate(x, y, gate);
                if (isInside && !gateNames.includes(gate.name)) {
                    newGateNames = [...gateNames, gate.name];
                    this.setGateNamesForCell(i, newGateNames);
                }
            } else {
                if (gateNames.includes(gate.name)) {
                    newGateNames = gateNames.filter((gateName) => gateName !== gate.name);
                    this.setGateNamesForCell(i, newGateNames);
                }
            }
        }
    }

    /**
     * Add a new gate to gates array and update all cells in the selection with gate value
     * Update the dataStore with the updated gates
     *
     * @param gate - Gate to be added
     */
    async addGate(gate: Gate) {
        // Ensure column data is loaded
        await this.ensureGateColumnDataLoaded();

        const clonedGate = {
            ...gate,
            geometry: JSON.parse(JSON.stringify(gate.geometry)), // Only clone geometry
        };
        this.gates.set(gate.id, clonedGate);
        this.updateCellsWithGate(clonedGate, true);
        this.markDirty();

        if (this.gateColumn) {
            this.dataStore.dataChanged([GATES_COLUMN_NAME]);
        }
    }

    /**
     * Getter for the mobx computed gatesArray
     * Deep clone the gates to avoid mutating mobx properties
     */
    get gatesArray(): Gate[] {
        return Array.from(this.gates.values()).map((gate) => ({
            ...gate,
            geometry: JSON.parse(JSON.stringify(gate.geometry)),
        }));
    }

    /**
     * Get gates for the columns
     * @param xField - X Column field of chart
     * @param yField - Y Column field of chart
     * @returns - Gates present in the chart with the given x and y param fields
     */
    getGatesForColumns(xField: string, yField: string): Gate[] {
        return this.gatesArray.filter((gate) => gate.columns[0] === xField && gate.columns[1] === yField);
    }

    /**
     * Check if gatesArray has the given gate name
     */
    hasGateName(name: string): boolean {
        return this.gatesArray.some((gate) => gate.name === name);
    }

    /**
     * Update the gates column when it's loaded
     * Ensure that the gates column and it's data is loaded before updating the column
     */
    public async updateGatesColumnWhenReady() {
        if (this.gates.size > 0) {
            // Ensure column data is loaded first
            await this.ensureGateColumnDataLoaded();
            await this.updateGateColumn();
        }
    }

    private splitGateNames(value?: string): string[] {
        if (!value) return [];
        return value.split(/\s*;\s*/).filter((v) => v.trim());
    }

    private joinGateNames(names: string[]): string {
        return names.filter((n) => n).join("; ");
    }

    /**
     * Updating the dataStore.config.gates and dataStore.gates with the updated gates
     * Storing gates in dirtyMetadata for persistence
     * Call the dataChanged function of dirtyColumns
     */
    private markDirty() {
        const gatesArray = this.gatesArray;
        this.dataStore.config.gates = gatesArray;
        this.dataStore.gates = gatesArray;
        this.dataStore.dirtyMetadata.add("gates");

        if (this.gateColumn) {
            (this.dataStore.dirtyColumns as any).data_changed[GATES_COLUMN_NAME] = true;
        }
    }

    toJSON(): Gate[] {
        return this.gatesArray;
    }
}
