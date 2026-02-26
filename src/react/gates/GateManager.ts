import type DataStore from "@/datastore/DataStore";
import type { Gate } from "./types";
import { action, computed, makeObservable, observable } from "mobx";
import type { LoadedDataColumn } from "@/charts/charts";
import { extractCoords, isPointInGate } from "./gateUtils";
import { loadColumn } from "@/dataloaders/DataLoaderUtil";

const GATES_COLUMN_NAME = "__gates__";
const GATE_NONE_VALUE = "N/A";
/**
 * Class which manages the gating operations
 * Creates a gates column if it doesn't exist, loads the gates column if it's not loaded
 * Contains the gate lifecycle methods for adding, updating and deleting
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
            updateGate: action,
            deleteGate: action,
            gatesArray: computed,
        });

        this.initializeGateColumn();
    }

    /**
     * Load the gates column if not loaded yet
     * Add the new gate to the gates map
     * Update all the cells inside the gate geometry with the gate name
     * Update the data store with the new gate
     * @param gate - new gate to be added
     */
    async addGate(gate: Gate) {
        // Ensure column data is loaded
        await this.ensureGatesColumnDataLoaded();

        const clonedGate = {
            ...gate,
            geometry: JSON.parse(JSON.stringify(gate.geometry)), // Only clone geometry
        };

        action(() => {
            this.gates.set(gate.id, clonedGate);
        })();

        // Update the cells with the gate
        this.updateCellsWithGate(clonedGate, true);
        this.updateDataStoreWithGates();

        if (this.gateColumn) {
            this.dataStore.dataChanged([GATES_COLUMN_NAME]);
        }
    }

    /**
     * Update the gate object with the updated gate
     * Update the cells inside the gate only if the geometry or name or columns change exist
     * Update the datastore
     * @param gateId - gate id
     * @param updates - updated gate object (partial or full)
     */
    updateGate(gateId: string, updates: Partial<Gate>) {
        const gate = this.gates.get(gateId);
        if (!gate) {
            console.error(`Gate ${gateId} not found`);
            return;
        }

        const updatedGate = {
            ...gate,
            ...updates,
            geometry: updates.geometry ? JSON.parse(JSON.stringify(updates.geometry)) : gate.geometry,
        };

        action(() => {
            this.gates.set(gateId, updatedGate);
        })();

        const cellUpdatedNeeded = "geometry" in updates || "name" in updates || "columns" in updates;

        if (cellUpdatedNeeded) {
            this.updateCellsWithGate(updatedGate, true);
        }

        this.updateDataStoreWithGates();

        if (this.gateColumn) {
            this.dataStore.dataChanged([GATES_COLUMN_NAME]);
        }
    }

    /**
     * Remove the gate name from the cells in gate column
     * Remove the gate from the gates object
     * Rebuild the values array to remove deleted values (required for filtering)
     * Update the datastore
     */
    deleteGate(gateId: string) {
        const gate = this.gates.get(gateId);
        if (!gate) {
            console.error(`Gate ${gateId} not found`);
            return;
        }

        this.updateCellsWithGate(gate, false);

        action(() => {
            this.gates.delete(gateId);
        })();

        this.rebuildValuesArray();
        this.updateDataStoreWithGates();

        if (this.gateColumn) {
            this.dataStore.dataChanged([GATES_COLUMN_NAME]);
        }
    }

    /**
     * Get or create the gate column
     * Initialise the gate column with empty value initially
     */
    private initializeGateColumn() {
        const existingCol = this.dataStore.columnIndex[GATES_COLUMN_NAME];
        if (existingCol && existingCol.datatype === "multitext") {
            // Get the existing column
            this.gateColumn = existingCol as LoadedDataColumn<"multitext">;
        } else {
            // Create a new column
            const column = {
                name: "gates",
                field: GATES_COLUMN_NAME,
                datatype: "multitext" as const,
                editable: false, // as of now
                delimiter: "," as const,
                values: [GATE_NONE_VALUE] as string[],
                stringLength: 10,
            };

            const data = new SharedArrayBuffer(this.dataStore.size * column.stringLength * 2);
            const dataArray = new Uint16Array(data);

            dataArray.fill(65535);

            // Mark the first value of all cells as 'N/A' initially
            for (let i=0; i<this.dataStore.size; i++) {
                dataArray[i * column.stringLength] = 0;
            }

            // Add the column to datastore
            this.dataStore.addColumn(column, data, true);

            this.gateColumn = this.dataStore.columnIndex[GATES_COLUMN_NAME] as LoadedDataColumn<"multitext">;
        }

        this.loadGatesFromConfig();
    }

    /**
     * Load the gates stored in the dataStore.config or dataStore
     * Rebuild the gates column with the stored gates
     */
    private async loadGatesFromConfig() {
        // Get the stored gates
        const savedGates = this.dataStore.config.gates || this.dataStore.gates;

        if (Array.isArray(savedGates)) {
            for (const gate of savedGates) {
                // Deep clone the gate
                const clonedGate = {
                    ...gate,
                    geometry: JSON.parse(JSON.stringify(gate.geometry)),
                };
                this.gates.set(gate.id, clonedGate);
            }
        }

        await this.rebuildGatesColumnWhenReady();
    }

    /**
     * Load the gates column if not loaded
     * Rebuild the gates column
     */
    public async rebuildGatesColumnWhenReady() {
        if (this.gates.size > 0) {
            // Ensure column data is loaded before updating gate column
            await this.ensureGatesColumnDataLoaded();
            await this.rebuildGateColumn();
        }
    }

    /**
     * Check if the gates column is loaded or not
     * Load the column explicitly if not loaded yet
     */
    private async ensureGatesColumnDataLoaded() {
        if (!this.gateColumn) {
            this.gateColumn = this.dataStore.columnIndex[GATES_COLUMN_NAME] as LoadedDataColumn<"multitext">;
            if (!this.gateColumn) {
                console.error("No gates column in the dataStore");
                return;
            }
        }

        if (this.gateColumn.data) {
            // Gate column data already loaded
            return;
        }

        // Explicitly load the column data
        try {
            await loadColumn(this.dataStore.name, GATES_COLUMN_NAME);
            this.gateColumn = this.dataStore.columnIndex[GATES_COLUMN_NAME] as LoadedDataColumn<"multitext">;
            if (!this.gateColumn.data) {
                console.error("Failed to load gates column data");
            }
        } catch (error) {
            console.error("Failed to load gates column data:", error);
        }
    }

    /**
     * Ensure the gate column data is loaded
     * Set the gate names as empty for all cells
     * Recompute gates for the cells
     */
    private async rebuildGateColumn() {
        if (!this.gateColumn) {
            return;
        }

        // Populate the gate names for all cells with empty array
        for (let i = 0; i < this.dataStore.size; i++) {
            this.setGateNamesForCell(i, [GATE_NONE_VALUE]);
        }

        // Recompute gate membership
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

        // Only notify if something was updated
        if (anyColumnsLoaded) {
            this.dataStore.dataChanged([GATES_COLUMN_NAME]);
        }
    }

    /**
     * Get gate names for each cell
     * To add the gate name to the cell, check if it inside and add the gate name
     * To remove, remove the gate from the cell if the cell already has the gate name
     * @param gate - Gate object
     * @param add - true if add, false if remove
     */
    private updateCellsWithGate(gate: Gate, add: boolean) {
        if (!this.gateColumn) return;

        const [xField, yField] = gate.columns;

        const xCol = this.dataStore.columnIndex[xField];
        const yCol = this.dataStore.columnIndex[yField];

        if (!xCol || !yCol || !xCol?.data || !yCol?.data) return;

        // Extract polygon coords from gate geometry
        const polygonCoords = extractCoords(gate.geometry);
        if (polygonCoords.length === 0) return;

        for (let i = 0; i < this.dataStore.size; i++) {
            const x = xCol.data[i];
            const y = yCol.data[i];

            if (!Number.isFinite(x) || !Number.isFinite(y)) continue;

            // Get gate names for the cell index - i
            const gateNames = this.getGateNamesForCell(i);

            let newGateNames: string[];

            if (add) {
                // Add the gate name to the cell if it is inside the gate
                const isInside = isPointInGate(x, y, gate);
                if (isInside && !gateNames.includes(gate.name)) {
                    const currentGates = gateNames.filter((g) => g !== GATE_NONE_VALUE);
                    newGateNames = [...currentGates, gate.name];
                    this.setGateNamesForCell(i, newGateNames);
                }
            } else {
                // Remove the gate name from the cell if it previously was inside the gate
                if (gateNames.includes(gate.name)) {
                    newGateNames = gateNames.filter((gateName) => gateName !== gate.name);
                    if (newGateNames.length === 0) newGateNames = [GATE_NONE_VALUE];
                    this.setGateNamesForCell(i, newGateNames);
                }
            }
        }
    }

    /**
     * Get the value index of the gates associated with the cells
     * @param cellIndex - Index of the cell in dataStore
     * @returns gate names for the cell
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
     * Assign the value indices of all gate names to the cell
     * @param cellIndex - Index of the cell in dataStore
     * @param gateNames - Gate names to update the cell with
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
     * Update the datastore and datastore.config with the gates
     * Add the gates column to dirtyMetaData to persist
     * Notify dataStore.dirtyColumns that data of this column is changed
     */
    public updateDataStoreWithGates() {
        const gatesArray = this.gatesArray;
        // Update dataStore and config with the new gates array
        this.dataStore.config.gates = gatesArray;
        this.dataStore.gates = gatesArray;
        this.dataStore.dirtyMetadata.add("gates");

        if (this.gateColumn) {
            (this.dataStore.dirtyColumns as any).data_changed[GATES_COLUMN_NAME] = true;
        }
    }

    /**
     * Get or add the index of the value of the gate name
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
     * Remove the gate name from the values array
     * This is required for the selection dialog to be in sync with the current gates
     * and not show the deleted gates
     * Reassign the indices of the cells with the updated values array so the cells 
     * have the latest index value of the updated values array
     */
    private rebuildValuesArray() {
        if (!this.gateColumn) return;
    
        // Read current gate names per cell
        const gateNamesPerCell: string[][] = [];
        for (let i = 0; i < this.dataStore.size; i++) {
            gateNamesPerCell[i] = this.getGateNamesForCell(i);
        }
    
        // Replace values with only current gate names (so selection dialog dropdown stays in sync)
        const currentGateNames = [GATE_NONE_VALUE, ...Array.from(this.gates.values()).map((g) => g.name).sort()];
        this.gateColumn.values.length = 0;
        this.gateColumn.values.push(...currentGateNames);
    
        // Write back each cell so indices point into the new values
        for (let i = 0; i < this.dataStore.size; i++) {
            this.setGateNamesForCell(i, gateNamesPerCell[i]);
        }
    }

    toJSON(): Gate[] {
        return this.gatesArray;
    }

    // Getter for the computed gatesArray, create a deep clone of geometry
    get gatesArray(): Gate[] {
        return Array.from(this.gates.values()).map((gate) => ({
            ...gate,
            geometry: JSON.parse(JSON.stringify(gate.geometry)),
        }));
    }

    getGatesForColumns(xField: string, yField: string): Gate[] {
        return this.gatesArray.filter((gate) => gate.columns[0] === xField && gate.columns[1] === yField);
    }

    hasGateName(name: string): boolean {
        return this.gatesArray.some((gate) => gate.name === name);
    }
}
