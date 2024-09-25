import type BaseChart from "@/charts/BaseChart";
import type { BaseDialog } from "@/utilities/Dialog";
import { observer } from "mobx-react-lite";
import { createContext, useContext, useEffect, useState } from "react";

/** All charts and dialogs have this as an `observable` property that can be passed to
 * {@link react_utils.createMdvPortal} to allow access to the container they should render things like
 * which allows {@link mobx-react-lite.observer} components to {@link useOuterContainer} for access to the container they should render things like
 * tooltips and dropdowns into.
 */
export type OuterContainerObservable = { container: HTMLElement };
const OuterContainerContext = createContext<Element>(null);
/**
 * As of now we only ever use react to render charts & dialogs, and they should be mobx observable
 * in as much as any `changeBaseDocument()` or `setParent()` will cause this to update...
 * hopefully the implementation will be easy to change if/when necessary.
 */
export const OuterContainerProvider = observer(
    ({
        children,
        parent,
    }: { children: JSX.Element; parent?: BaseChart | BaseDialog }) => {
        return (
            <OuterContainerContext.Provider
                value={parent?.observable?.container || document.body}
            >
                {children}
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
