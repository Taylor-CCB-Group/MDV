import { useEffect, useLayoutEffect, useRef, useState } from "react";
import { renderMermaidDiagram } from "./mermaidUtils";
import { registerTextBoxFencedRenderer } from "./textBoxFencedRendererRegistry";

type MermaidState = {
    error: string | null;
    svg: string | null;
};

export function MermaidDiagramBlock({ code }: { code: string }) {
    const blockRef = useRef<HTMLDivElement | null>(null);
    const [state, setState] = useState<MermaidState>({ svg: null, error: null });
    const [htmlClassRevision, setHtmlClassRevision] = useState(0);

    useLayoutEffect(() => {
        const el = blockRef.current;
        if (!el) {
            return;
        }
        const root = el.ownerDocument.documentElement;
        const obs = new MutationObserver(() => {
            setHtmlClassRevision((n) => n + 1);
        });
        obs.observe(root, { attributes: true, attributeFilter: ["class"] });
        return () => {
            obs.disconnect();
        };
    }, []);

    useEffect(() => {
        let cancelled = false;

        setState({ svg: null, error: null });
        const targetDoc = blockRef.current?.ownerDocument ?? document;
        renderMermaidDiagram(code, targetDoc)
            .then((svg) => {
                if (cancelled) return;
                setState({ svg, error: null });
            })
            .catch((error) => {
                if (cancelled) return;
                const message = error instanceof Error ? error.message : String(error);
                setState({
                    svg: null,
                    error: message,
                });
            });

        return () => {
            cancelled = true;
        };
    }, [code, htmlClassRevision]);

    return (
        <div className="mdv-textbox-mermaid-block" ref={blockRef}>
            {state.svg ? (
                <div
                    className="mdv-textbox-mermaid-diagram"
                    // biome-ignore lint/security/noDangerouslySetInnerHtml: sanitized in mermaidUtils before rendering
                    dangerouslySetInnerHTML={{ __html: state.svg }}
                />
            ) : null}
            {state.error ? (
                <>
                    <div className="mdv-textbox-mermaid-error">
                        Unable to render Mermaid diagram.
                    </div>
                    <pre className="mdv-textbox-mermaid-fallback">
                        <code>{code}</code>
                    </pre>
                </>
            ) : null}
            {!state.svg && !state.error ? (
                <div className="mdv-textbox-mermaid-loading">Rendering diagram…</div>
            ) : null}
        </div>
    );
}

registerTextBoxFencedRenderer("mermaid", MermaidDiagramBlock);
