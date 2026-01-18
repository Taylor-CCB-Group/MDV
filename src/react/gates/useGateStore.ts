import { useMemo } from 'react';
import { useDataStore } from '../context';
import { getGateStore } from './gateStoreManager';

/**
 * Hook to get the GateStore for the current chart's datasource
 */
export function useGateStore() {
    const dataStore = useDataStore();
    return useMemo(() => getGateStore(dataStore), [dataStore]);
}

