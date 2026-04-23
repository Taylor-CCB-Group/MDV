let mermaidIdCounter = 0;

type ThemeMode = "dark" | "light";

function getCurrentThemeMode(doc: Document): ThemeMode {
    return doc.documentElement.classList.contains("dark") ? "dark" : "light";
}

function getThemeConfig(doc: Document) {
    const rootStyle = getComputedStyle(doc.documentElement);
    const themeMode = getCurrentThemeMode(doc);
    const textColor = rootStyle.getPropertyValue("--text_color").trim() || (themeMode === "dark" ? "#ffffff" : "#111111");
    const mainPanelColor =
        rootStyle.getPropertyValue("--main_panel_color").trim() ||
        (themeMode === "dark" ? "#000000" : "#f1f1f1");
    const primaryBackground =
        rootStyle.getPropertyValue("--primary_background").trim() ||
        (themeMode === "dark" ? "#353536" : "#ecf0f1");
    const borderColor =
        rootStyle.getPropertyValue("--border_menu_bar_color").trim() ||
        (themeMode === "dark" ? "#727272" : "#c9c9c9");

    return {
        startOnLoad: false,
        securityLevel: "strict" as const,
        theme: "base" as const,
        darkMode: themeMode === "dark",
        themeVariables: {
            background: "transparent",
            fontFamily: "inherit",
            primaryColor: mainPanelColor,
            primaryTextColor: textColor,
            primaryBorderColor: borderColor,
            lineColor: textColor,
            textColor,
            mainBkg: mainPanelColor,
            secondBkg: primaryBackground,
            tertiaryColor: primaryBackground,
            clusterBkg: primaryBackground,
            clusterBorder: borderColor,
            nodeBorder: borderColor,
            edgeLabelBackground: mainPanelColor,
        },
        themeCSS: `
            svg { background-color: transparent !important; }
            .label text, .nodeLabel text, .edgeLabel text, .cluster-label text, .flowchart-label text, .label text tspan {
                fill: ${textColor} !important;
            }
            .label, .nodeLabel, .edgeLabel, .cluster-label, .flowchart-label {
                color: ${textColor} !important;
            }
            .label foreignObject div, .label foreignObject span, .nodeLabel p, .nodeLabel span, .nodeLabel div {
                color: ${textColor} !important;
            }
            .edgePath path, .flowchart-link {
                stroke: ${textColor} !important;
                fill: none !important;
            }
            .arrowheadPath, .marker path, marker path {
                stroke: ${textColor} !important;
                fill: ${textColor} !important;
            }
            .node rect, .node circle, .node ellipse, .node polygon, .node path {
                fill: ${primaryBackground} !important;
                stroke: ${borderColor} !important;
            }
            .edgeLabel rect, .labelBkg {
                fill: ${mainPanelColor} !important;
                opacity: 0.95 !important;
            }
            .cluster rect {
                fill: ${primaryBackground} !important;
                stroke: ${borderColor} !important;
            }
        `,
    };
}

async function getMermaidRuntime() {
    const mermaidModule = await import("mermaid");
    return mermaidModule.default;
}

function stripEventHandlerAttributes(root: Element) {
    const walk = (el: Element) => {
        for (const { name } of [...el.attributes]) {
            if (/^on/i.test(name)) {
                el.removeAttribute(name);
            }
        }
        for (const child of el.children) {
            walk(child);
        }
    };
    walk(root);
}

/**
 * Mermaid node labels are HTML inside <foreignObject>. A single SVG DOMPurify pass can
 * leave those tags in place but drop their inner content. Sanitize XHTML in each
 * foreignObject with the HTML profile, then harden the tree (no scripts, no on*).
 */
async function sanitizeMermaidSvgString(svg: string) {
    const domPurifyModule = await import("dompurify");
    const purify = domPurifyModule.default;
    const parser = new DOMParser();
    const parsed = parser.parseFromString(svg, "image/svg+xml");
    if (parsed.querySelector("parsererror") !== null) {
        return purify.sanitize(svg, {
            USE_PROFILES: { svg: true, svgFilters: true },
        });
    }
    const root = parsed.documentElement;
    for (const node of root.querySelectorAll("script")) {
        node.remove();
    }
    for (const fo of root.querySelectorAll("foreignObject")) {
        const inner = fo.innerHTML;
        if (inner) {
            fo.innerHTML = purify.sanitize(inner, { USE_PROFILES: { html: true } });
        }
    }
    stripEventHandlerAttributes(root);
    return new XMLSerializer().serializeToString(root);
}

export async function renderMermaidDiagram(code: string, doc: Document = document) {
    const mermaid = await getMermaidRuntime();
    const theme = getThemeConfig(doc);
    mermaid.initialize(theme);
    const diagramId = `mdv-textbox-mermaid-${mermaidIdCounter++}`;
    const result = await mermaid.render(diagramId, code);
    return sanitizeMermaidSvgString(result.svg);
}

export function resetMermaidRuntimeForTests() {
    mermaidIdCounter = 0;
}
