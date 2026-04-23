import { createElement, useEffect, useMemo, useState, type ReactNode } from "react";
import remarkCollapsibleHeadings from "./remarkCollapsibleHeadings";
import { getTextBoxFencedRenderer } from "./textBoxFencedRendererRegistry";
import "./MermaidDiagramBlock";

type MarkdownRuntime = {
    ReactMarkdown: any;
    remarkGfm: any;
    rehypeRaw: any;
    rehypeSanitize: any;
};

function getCodeBlockLanguage(className: string | undefined) {
    if (!className) return "";
    const languageMatch = /language-([\w-]+)/.exec(className);
    return languageMatch ? languageMatch[1].toLowerCase() : "";
}

function toCodeString(children: ReactNode) {
    if (Array.isArray(children)) return children.join("");
    return String(children ?? "");
}

const loadingMarkup = "Rendering markdown…";

const TextBoxMarkdownRenderer = ({ markdown }: { markdown: string }) => {
    const [runtime, setRuntime] = useState<MarkdownRuntime | null>(null);

    useEffect(() => {
        let cancelled = false;
        Promise.all([
            import("react-markdown"),
            import("remark-gfm"),
            import("rehype-raw"),
            import("rehype-sanitize"),
        ])
            .then(
                ([
                    reactMarkdownModule,
                    remarkGfmModule,
                    rehypeRawModule,
                    rehypeSanitizeModule,
                ]) => {
                    if (cancelled) return;
                    setRuntime({
                        ReactMarkdown: reactMarkdownModule.default,
                        remarkGfm: remarkGfmModule.default,
                        rehypeRaw: rehypeRawModule.default,
                        rehypeSanitize: rehypeSanitizeModule.default,
                    });
                },
            )
            .catch((error) => {
                if (cancelled) return;
                console.error("failed to load markdown runtime", error);
            });

        return () => {
            cancelled = true;
        };
    }, []);

    const markdownComponents = useMemo(
        () => ({
            pre({
                children,
                className,
            }: {
                children?: ReactNode;
                className?: string;
            }) {
                return (
                    <div
                        className={
                            className
                                ? `mdv-textbox-fenced ${className}`
                                : "mdv-textbox-fenced"
                        }
                    >
                        {children}
                    </div>
                );
            },
            code({
                className,
                children,
                node: _node,
                ...props
            }: {
                className?: string;
                children?: ReactNode;
                node?: unknown;
            }) {
                const language = getCodeBlockLanguage(className);
                const renderer = language
                    ? getTextBoxFencedRenderer(language)
                    : undefined;
                if (!renderer) {
                    const blockLang = getCodeBlockLanguage(className);
                    const text = toCodeString(children).replace(/\n$/, "");
                    const isMultiline = text.includes("\n");
                    if (!isMultiline) {
                        return (
                            <code className={className} {...props}>
                                {children}
                            </code>
                        );
                    }
                    const summaryLabel = blockLang || "text";
                    return (
                        <details className="mdv-textbox-codeblock" open>
                            <summary className="mdv-textbox-summary mdv-textbox-codeblock-summary">
                                {summaryLabel}
                            </summary>
                            <code className={className} {...props}>
                                {children}
                            </code>
                        </details>
                    );
                }

                const code = toCodeString(children).replace(/\n$/, "");
                return createElement(renderer, { code });
            },
        }),
        [],
    );

    if (!runtime) {
        return (
            <div className="markdown-body mdv-textbox-markdown">
                <div className="mdv-textbox-loading">{loadingMarkup}</div>
            </div>
        );
    }

    const { ReactMarkdown, remarkGfm, rehypeRaw, rehypeSanitize } = runtime;

    return (
        <div className="markdown-body mdv-textbox-markdown">
            <ReactMarkdown
                remarkPlugins={[remarkGfm, remarkCollapsibleHeadings]}
                rehypePlugins={[rehypeRaw, rehypeSanitize]}
                components={markdownComponents}
            >
                {markdown}
            </ReactMarkdown>
        </div>
    );
};

export default TextBoxMarkdownRenderer;
