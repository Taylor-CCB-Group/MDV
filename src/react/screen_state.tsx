import type BaseChart from "@/charts/BaseChart";
import type { BaseDialog } from "@/utilities/Dialog";
import { observer } from "mobx-react-lite";
import { createContext, useContext, type PropsWithChildren } from "react";
import type { createMdvPortal } from "./react_utils";
import createCache, { type EmotionCache } from "@emotion/cache";
import { CacheProvider } from "@emotion/react";

/** All charts and dialogs have this as an `observable` property that can be passed to
 * {@link createMdvPortal} to allow access to the container they should render things like
 * which allows Mobx Observer components to {@link useOuterContainer} for access to the container they should render things like
 * tooltips and dropdowns into.
 */
export type OuterContainerObservable = { container: HTMLElement };
const OuterContainerContext = createContext<HTMLElement>(null as any);

let idCounter = 0;

function numberToAlpha(n: number): string {
    let result = "";
    while (n > 0) {
      const remainder = (n - 1) % 26;
      result = String.fromCharCode(97 + remainder) + result;
      n = Math.floor((n - 1) / 26);
    }
    return result;
}

// unique key generator
export function generateCacheKey() {
  idCounter++;
  return `mui-${numberToAlpha(idCounter)}`;
}

// Cache map for storing EmotionCache for separate windows
const cacheMap = new Map<HTMLElement, EmotionCache>();

// Getting or creating EmotionCache for storing cache using CacheProvider
export function useEmotionCache(container: HTMLElement) {
    let cache = cacheMap.get(container);
    if (!cache) {
        cache = createCache({ key: generateCacheKey(), container });
        cacheMap.set(container, cache);
    }
    return cache;
}
/**
 * As of now we only ever use react to render charts & dialogs, and they should be mobx observable
 * in as much as any `changeBaseDocument()` or `setParent()` will cause this to update...
 * hopefully the implementation will be easy to change if/when necessary.
 */
export const OuterContainerProvider = observer(
    ({ children, parent }: { parent?: BaseChart<any> | BaseDialog } & PropsWithChildren ) => {
        const container = parent?.observable?.container || document.body;
        const cache = useEmotionCache(container);
        return (
            <OuterContainerContext.Provider value={container}>
                <CacheProvider value={cache}>
                    {children}
                </CacheProvider>
            </OuterContainerContext.Provider>
        );
    },
);
/**
 * This should allow any component that may need knowledge of an outer container to render into
 * (e.g. for displaying a tooltip or dropdown that may not fit inside that part of DOM).
 *
 * This should consistently provide the right result with any permutation of fullscreen & popout.
 *
 * Specification of how this is implemented subject to change, but it should *always* be valid to call this
 * from any MDV react component (there may be situations in which it doesn't return anything more useful than `document.body`).
 */
export const useOuterContainer = () => {
    const container = useContext(OuterContainerContext);
    return container || document.body; //this wasn't getting the right default when leaving fullscreen (d.setParent(null))
};
