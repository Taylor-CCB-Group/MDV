import type DataStore from "@/datastore/DataStore";
import type { Gate } from "./types";
import { action, computed, makeObservable, observable } from "mobx";
import type { LoadedDataColumn } from "@/charts/charts";
import { extractCoords, isPointInGate } from "./gateUtils";
import { loadColumn } from "@/dataloaders/DataLoaderUtil";

// todo: Change later
const GATES_COLUMN_NAME = "__gates__";
export class GateStore {
    gates: Map<string, Gate> = new Map();

    private dataStore: DataStore;
    private gateColumn: LoadedDataColumn<"multitext"> | null = null;

    constructor(dataStore: DataStore) {
        this.dataStore = dataStore;

        makeObservable(this, {
            gates: observable,
            addGate: action,
            // removeGate: action,
            // updateGate: action,
            gatesArray: computed,
        });

        this.initializeGateColumn();
        this.loadGatesFromConfig();
    }

    private async ensureGateColumnDataLoaded() {
        if (!this.gateColumn) return;

        // Check if data is already loaded
        if (this.gateColumn.data) {
            return; // Already loaded
        }

        // Check if column exists in columnIndex
        if (!this.dataStore.columnIndex[GATES_COLUMN_NAME]) {
            return; // Column doesn't exist yet
        }

        // Explicitly load the column data
        try {
            await loadColumn(this.dataStore.name, GATES_COLUMN_NAME);
            // After loading, the gateColumn reference should now have data
            this.gateColumn = this.dataStore.columnIndex[GATES_COLUMN_NAME] as LoadedDataColumn<"multitext">;
        } catch (error) {
            console.error("Failed to load gates column data:", error);
        }
    }

    private initializeGateColumn() {
        const existingCol = this.dataStore.columnIndex[GATES_COLUMN_NAME];
        if (existingCol && existingCol.datatype === "multitext") {
            this.gateColumn = existingCol as LoadedDataColumn<"multitext">;

            // CRITICAL: Ensure the column has data loaded
            if (!this.gateColumn.data) {
                // Column exists but data isn't loaded - this shouldn't happen for multitext
                // but if it does, we need to wait for it to load
                console.warn("Gates column exists but data is not loaded yet");
            }

            // Ensure stringLength is set if missing
            if (typeof this.gateColumn.stringLength !== "number") {
                // @ts-expect-error: Forcibly assign stringLength if not set
                this.gateColumn.stringLength = 10 as any;
            }

            // Ensure values array exists
            if (!Array.isArray(this.gateColumn.values)) {
                this.gateColumn.values = [];
            }

            // If data isn't loaded, we'll load it asynchronously
            // Don't block here, but ensure it gets loaded
            if (!this.gateColumn.data) {
                // Load it asynchronously - don't await, but ensure it happens
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

    private async loadGatesFromConfig() {
        const savedGates = this.dataStore.config.gates || this.dataStore.gates;
        if (Array.isArray(savedGates)) {
            for (const gate of savedGates) {
                // this.gates.set(gate.id, gate);
                const clonedGate = {
                    ...gate,
                    geometry: JSON.parse(JSON.stringify(gate.geometry)),
                };
                this.gates.set(gate.id, clonedGate);
            }
        }
        // Ensure column data is loaded before updating
        await this.ensureGateColumnDataLoaded();

        // Now try to update (will only work if required columns are loaded)
        this.updateGatesColumnWhenReady().catch(console.error);
    }

    private getValueIndex(gateName: string): number {
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

    private setGateNamesForCell(cellIndex: number, gateNames: string[]) {
        if (!this.gateColumn || !this.gateColumn.stringLength) return;

        const baseIndex = cellIndex * this.gateColumn.stringLength;
        const maxValues = Math.min(gateNames.length, this.gateColumn.stringLength);

        for (let i = 0; i < this.gateColumn.stringLength; i++) {
            this.gateColumn.data[baseIndex + i] = 65535;
        }

        const sortedNames = [...gateNames].sort();

        for (let i = 0; i < maxValues; i++) {
            const valueIndex = this.getValueIndex(sortedNames[i]);
            this.gateColumn.data[baseIndex + i] = valueIndex;
        }
    }

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
            // If columns aren't loaded yet, we'll update later when they are
        }

        // Only notify if we actually updated something
        if (anyColumnsLoaded) {
            this.dataStore.dataChanged([GATES_COLUMN_NAME]);
        }
    }

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

    // removeGate(gateId: string) {
    //     const gate = this.gates.get(gateId);
    //     if (gate) {
    //         this.updateCellsWithGate(gate, false);
    //     }

    //     this.gates.delete(gateId);
    //     this.markDirty();

    //     if (this.gateColumn) this.dataStore.dataChanged([GATES_COLUMN_NAME]);
    // }

    // updateGate(gateId: string, updates: Partial<Gate>) {
    //     const gate = this.gates.get(gateId);
    //     if (gate) {
    //         this.gates.set(gateId, { ...gate, ...updates });
    //         this.markDirty();
    //     }
    // }

    get gatesArray(): Gate[] {
        // return Array.from(this.gates.values());
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
