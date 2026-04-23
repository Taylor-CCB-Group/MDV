import { useEffect, useState } from "react";
import { renderMermaidDiagram } from "./mermaidUtils";
import { registerTextBoxFencedRenderer } from "./textBoxFencedRendererRegistry";

type MermaidState = {
    error: string | null;
    svg: string | null;
};

export function MermaidDiagramBlock({ code }: { code: string }) {
    const [state, setState] = useState<MermaidState>({ svg: null, error: null });

    useEffect(() => {
        let cancelled = false;

        setState({ svg: null, error: null });
        renderMermaidDiagram(code)
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
    }, [code]);

    return (
        <div className="mdv-textbox-mermaid-block">
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
