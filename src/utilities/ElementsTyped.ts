import { addElProps } from "./Elements";

export type TagKey = keyof HTMLElementTagNameMap;

export type Attrs = {
    styles?: Partial<CSSStyleDeclaration>;
    classes?: string[];
    text?: string;
    [key: string]: string | string[] | Partial<CSSStyleDeclaration> | undefined;
};

/** Create an HTML element of the given type */
export function createEl<T extends TagKey>(
    type: T,
    attrs: Attrs,
    parent?: Element,
) {
    const el = document.createElement(type);

    if (attrs) {
        addElProps(el, attrs);
    }
    if (parent) {
        parent.append(el);
    }
    return el;
}
declare global {
    interface Element {
        draggable?: boolean;
    }
}