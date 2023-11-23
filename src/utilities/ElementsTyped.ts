import { addElProps } from "./Elements";

type TagKey = keyof HTMLElementTagNameMap

type Attrs = {
    styles?: Partial<CSSStyleDeclaration>,
    classes?: string[],
    text?: string,
    [key: string]: string | string[] | Partial<CSSStyleDeclaration> | undefined
}

/** Create an HTML element of the given type */
export function createEl<T extends TagKey>(type: T, attrs: Attrs, parent?: HTMLElement) {
    const el = document.createElement(type);

    if (attrs) {
        addElProps(el, attrs)
    }
    if (parent) {
        parent.append(el);
    }
    return el;
}