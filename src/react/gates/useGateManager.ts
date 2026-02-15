import { useMemo } from "react";
import { useDataStore } from "../context";
import type DataStore from "@/datastore/DataStore";
import { GateManager } from "./GateManager";

/**
 * Global map of DataStore -> GateManager
 * One GateManager per DataStore
 */
const gateManagers = new WeakMap<DataStore, GateManager>();

/**
 * Hook to get the GateManager for the current chart's datasource
 */
export function useGateManager() {
    const dataStore = useDataStore();

    return useMemo(() => {
        let gateManager = gateManagers.get(dataStore);

        if (!gateManager) {
            gateManager = new GateManager(dataStore);
            gateManagers.set(dataStore, gateManager);
        }

        return gateManager;
    }, [dataStore]);
}
