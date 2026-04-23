import { render, screen, waitFor } from "@testing-library/react";
import { beforeEach, describe, expect, test, vi } from "vitest";
import TextBoxMarkdownRenderer from "@/react/components/textbox/TextBoxMarkdownRenderer";
import { resetMermaidRuntimeForTests } from "@/react/components/textbox/mermaidUtils";

const mermaidMocks = vi.hoisted(() => ({
    initialize: vi.fn(),
    render: vi.fn(),
}));

const sanitizeMock = vi.hoisted(() => vi.fn((value: string) => value));

vi.mock("mermaid", () => ({
    default: {
        initialize: mermaidMocks.initialize,
        render: mermaidMocks.render,
    },
}));

vi.mock("dompurify", () => ({
    default: {
        sanitize: sanitizeMock,
    },
}));

describe("TextBoxMarkdownRenderer", () => {
    beforeEach(() => {
        mermaidMocks.initialize.mockReset();
        mermaidMocks.render.mockReset();
        sanitizeMock.mockClear();
        resetMermaidRuntimeForTests();
    });

    test("renders mermaid fenced blocks as diagrams", async () => {
        mermaidMocks.render.mockResolvedValue({
            svg: "<svg><text>diagram ok</text></svg>",
        });

        render(
            <TextBoxMarkdownRenderer
                markdown={"```mermaid\ngraph TD\nA-->B\n```"}
            />,
        );

        await waitFor(() => {
            expect(mermaidMocks.render).toHaveBeenCalled();
        });
        expect(mermaidMocks.initialize).toHaveBeenCalledWith(
            expect.objectContaining({
                theme: "base",
                themeCSS: expect.stringContaining("fill: none !important;"),
            }),
        );
        await waitFor(() => {
            const diagram = document.querySelector(".mdv-textbox-mermaid-diagram");
            expect(diagram?.innerHTML).toContain("<svg>");
        });
        expect(
            document.querySelector(".mdv-textbox-mermaid-fallback"),
        ).toBeNull();
    });

    test("falls back to raw code when mermaid rendering fails", async () => {
        mermaidMocks.render.mockRejectedValue(new Error("parse failure"));

        render(
            <TextBoxMarkdownRenderer
                markdown={"```mermaid\ngraph TD\nA-->B\n```"}
            />,
        );

        await waitFor(() => {
            expect(
                screen.getByText("Unable to render Mermaid diagram."),
            ).toBeTruthy();
        });
        const fallback = document.querySelector(".mdv-textbox-mermaid-fallback");
        expect(fallback?.textContent).toContain("graph TD");
        expect(fallback?.textContent).toContain("A-->B");
    });

    test("renders raw html content in markdown text", async () => {
        render(
            <TextBoxMarkdownRenderer
                markdown={"<h2>Legacy Header</h2><p><b>Legacy body</b></p>"}
            />,
        );

        await waitFor(() => {
            expect(screen.getByText("Legacy Header")).toBeTruthy();
        });
        expect(screen.getByText("Legacy body")).toBeTruthy();
    });
});
