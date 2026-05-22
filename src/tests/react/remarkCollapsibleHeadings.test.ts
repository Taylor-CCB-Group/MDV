import { describe, expect, test } from "vitest";
import { transformNodesForCollapsibleHeadings } from "@/react/components/textbox/remarkCollapsibleHeadings";

type RemarkNode = {
    type: string;
    depth?: number;
    value?: string;
    children?: RemarkNode[];
};

function heading(depth: number, text: string): RemarkNode {
    return {
        type: "heading",
        depth,
        children: [{ type: "text", value: text }],
    };
}

function paragraph(text: string): RemarkNode {
    return {
        type: "paragraph",
        children: [{ type: "text", value: text }],
    };
}

function matchingDetailsCloseIndex(nodes: RemarkNode[], openIndex: number) {
    let depth = 0;
    for (let index = openIndex; index < nodes.length; index += 1) {
        const node = nodes[index];
        if (node.type !== "html" || typeof node.value !== "string") continue;
        if (node.value.includes("<details")) depth += 1;
        if (node.value.includes("</details>")) {
            depth -= 1;
            if (depth === 0) return index;
        }
    }
    return -1;
}

describe("remarkCollapsibleHeadings", () => {
    test("wraps h2+ headings in open details with hierarchical boundaries", () => {
        const transformed = transformNodesForCollapsibleHeadings([
            heading(1, "Top"),
            paragraph("Intro"),
            heading(2, "Section A"),
            paragraph("A body"),
            heading(3, "Nested A"),
            paragraph("Nested body"),
            heading(2, "Section B"),
            paragraph("B body"),
        ]);

        expect(transformed[0]).toMatchObject({ type: "heading", depth: 1 });

        const sectionASummaryIndex = transformed.findIndex(
            (node) =>
                node.type === "html" &&
                typeof node.value === "string" &&
                node.value.includes("Section A"),
        );
        const sectionBSummaryIndex = transformed.findIndex(
            (node) =>
                node.type === "html" &&
                typeof node.value === "string" &&
                node.value.includes("Section B"),
        );
        const nestedSummaryIndex = transformed.findIndex(
            (node) =>
                node.type === "html" &&
                typeof node.value === "string" &&
                node.value.includes("Nested A"),
        );

        expect(sectionASummaryIndex).toBeGreaterThan(-1);
        expect(sectionBSummaryIndex).toBeGreaterThan(-1);
        expect(nestedSummaryIndex).toBeGreaterThan(-1);

        const sectionAOpenIndex = transformed.findIndex(
            (node, index) =>
                index <= sectionASummaryIndex &&
                node.type === "html" &&
                typeof node.value === "string" &&
                node.value.includes("<details open") &&
                node.value.includes("depth-2"),
        );
        const sectionACloseIndex = matchingDetailsCloseIndex(
            transformed,
            sectionAOpenIndex,
        );

        expect(sectionAOpenIndex).toBeGreaterThan(-1);
        expect(sectionACloseIndex).toBeGreaterThan(sectionAOpenIndex);
        expect(nestedSummaryIndex).toBeGreaterThan(sectionAOpenIndex);
        expect(nestedSummaryIndex).toBeLessThan(sectionACloseIndex);
        expect(sectionBSummaryIndex).toBeGreaterThan(sectionACloseIndex);
    });
});
