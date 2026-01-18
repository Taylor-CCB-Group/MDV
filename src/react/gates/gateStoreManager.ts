import { GateStore } from './GateStore';
import type DataStore from '@/datastore/DataStore';

/**
 * Global map of DataStore -> GateStore
 * One GateStore per DataStore
 */
const gateStores = new WeakMap<DataStore, GateStore>();

/**
 * Get or create a GateStore for a DataStore
 */
export function getGateStore(dataStore: DataStore): GateStore {
    let gateStore = gateStores.get(dataStore);
    
    if (!gateStore) {
        gateStore = new GateStore(dataStore);
        gateStores.set(dataStore, gateStore);
    }
    
    return gateStore;
}

