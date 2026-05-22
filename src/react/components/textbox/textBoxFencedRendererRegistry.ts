import type { ComponentType } from "react";

export type TextBoxFencedRendererProps = {
    code: string;
};

export type TextBoxFencedRenderer = ComponentType<TextBoxFencedRendererProps>;

const rendererRegistry = new Map<string, TextBoxFencedRenderer>();

export function registerTextBoxFencedRenderer(
    language: string,
    renderer: TextBoxFencedRenderer,
) {
    rendererRegistry.set(language.toLowerCase(), renderer);
}

export function getTextBoxFencedRenderer(language: string) {
    return rendererRegistry.get(language.toLowerCase());
}
