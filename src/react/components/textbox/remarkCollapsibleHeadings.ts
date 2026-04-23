type RemarkNode = {
    type: string;
    children?: RemarkNode[];
    depth?: number;
    value?: string;
};

function isHeadingNode(
    node: RemarkNode | undefined,
): node is RemarkNode & { depth: number; children: RemarkNode[] } {
    return Boolean(
        node &&
            node.type === "heading" &&
            typeof node.depth === "number" &&
            Array.isArray(node.children),
    );
}

function htmlNode(value: string): RemarkNode {
    return {
        type: "html",
        value,
    };
}

function escapeHtml(text: string) {
    return text
        .replaceAll("&", "&amp;")
        .replaceAll("<", "&lt;")
        .replaceAll(">", "&gt;")
        .replaceAll('"', "&quot;")
        .replaceAll("'", "&#39;");
}

function textFromNode(node: RemarkNode): string {
    if (typeof node.value === "string") return node.value;
    if (!Array.isArray(node.children)) return "";
    return node.children.map(textFromNode).join("");
}

function findSectionEnd(nodes: RemarkNode[], startIndex: number, depth: number) {
    let index = startIndex;
    while (index < nodes.length) {
        const candidate = nodes[index];
        if (isHeadingNode(candidate) && candidate.depth <= depth) {
            return index;
        }
        index += 1;
    }
    return nodes.length;
}

export function transformNodesForCollapsibleHeadings(nodes: RemarkNode[]) {
    const transformed: RemarkNode[] = [];
    let index = 0;

    while (index < nodes.length) {
        const node = nodes[index];
        if (!isHeadingNode(node) || node.depth < 2) {
            transformed.push(node);
            index += 1;
            continue;
        }

        const sectionEnd = findSectionEnd(nodes, index + 1, node.depth);
        const summaryText = textFromNode(node).trim() || `Section ${index + 1}`;
        const sectionChildren = transformNodesForCollapsibleHeadings(
            nodes.slice(index + 1, sectionEnd),
        );

        transformed.push(
            htmlNode(
                `<details open class="mdv-textbox-section mdv-textbox-section-depth-${node.depth}">`,
            ),
        );
        transformed.push(
            htmlNode(
                `<summary class="mdv-textbox-summary">${escapeHtml(summaryText)}</summary>`,
            ),
        );
        transformed.push(...sectionChildren);
        transformed.push(htmlNode("</details>"));

        index = sectionEnd;
    }

    return transformed;
}

export default function remarkCollapsibleHeadings() {
    return (tree: RemarkNode) => {
        if (!Array.isArray(tree.children)) return;
        tree.children = transformNodesForCollapsibleHeadings(tree.children);
    };
}
