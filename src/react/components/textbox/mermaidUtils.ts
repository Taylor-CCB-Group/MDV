let isMermaidInitialised = false;
let mermaidIdCounter = 0;

async function getMermaidRuntime() {
    const mermaidModule = await import("mermaid");
    const mermaid = mermaidModule.default;

    if (!isMermaidInitialised) {
        mermaid.initialize({
            startOnLoad: false,
            securityLevel: "strict",
        });
        isMermaidInitialised = true;
    }

    return mermaid;
}

async function sanitizeSvg(svg: string) {
    const domPurifyModule = await import("dompurify");
    const domPurify = domPurifyModule.default;
    return domPurify.sanitize(svg, {
        USE_PROFILES: {
            svg: true,
            svgFilters: true,
        },
    });
}

export async function renderMermaidDiagram(code: string) {
    const mermaid = await getMermaidRuntime();
    const diagramId = `mdv-textbox-mermaid-${mermaidIdCounter++}`;
    const result = await mermaid.render(diagramId, code);
    return sanitizeSvg(result.svg);
}

export function resetMermaidRuntimeForTests() {
    isMermaidInitialised = false;
    mermaidIdCounter = 0;
}
