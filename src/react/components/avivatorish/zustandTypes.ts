// can't see how to properly type things without copy/pasting from zustand
import type { StoreApi } from "zustand";

export type ExtractState<S> = S extends {
    getState: () => infer T;
}
    ? T
    : never;
export type ReadonlyStoreApi<T> = Pick<StoreApi<T>, "getState" | "subscribe">;
export type WithReact<S extends ReadonlyStoreApi<unknown>> = S & {
    getServerState?: () => ExtractState<S>;
};

// not copied from zustand source
export type ZustandStore<T> = WithReact<StoreApi<T>>;
export type Selector<S, U> = (state: ExtractState<S>) => U;
export type EqFn<U> = (a: U, b: U) => boolean;
