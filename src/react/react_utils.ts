import { createRoot } from "react-dom/client";

/**
 * todo - this is a placeholder for refactoring so that there is a single react root,
 * then charts/dialogs etc will be rendered into it with portals.
 */

const createMdvPortal = (component: JSX.Element, container: HTMLElement) => {
    const root = createRoot(container);
    root.render(component);
    return root;
}


export { createMdvPortal }