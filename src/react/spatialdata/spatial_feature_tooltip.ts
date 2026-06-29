import DOMPurify from "dompurify";
import type { SpatialFeatureTooltipData } from "@spatialdata/vis";

export function formatSpatialFeatureTooltipHtml(tooltip: SpatialFeatureTooltipData): string {
    const lines: string[] = [];
    if (tooltip.title) {
        lines.push(`<strong>${escapeHtml(tooltip.title)}</strong>`);
    }
    const sections = tooltip.sections ?? [];
    for (const section of sections) {
        if (section.title) {
            lines.push(`<div style="margin-top:6px;font-weight:600">${escapeHtml(section.title)}</div>`);
        }
        for (const item of section.items ?? []) {
            lines.push(
                `<div>${escapeHtml(item.label)}: ${escapeHtml(String(item.value ?? ""))}</div>`,
            );
        }
    }
    for (const item of tooltip.items ?? []) {
        lines.push(
            `<div>${escapeHtml(item.label)}: ${escapeHtml(String(item.value ?? ""))}</div>`,
        );
    }
    // Central sanitization of the assembled markup is the actual XSS guard for the
    // `innerHTML`-bound tooltip; `escapeHtml` below stays as defence in depth on each field.
    return DOMPurify.sanitize(lines.join(""));
}

function escapeHtml(value: string): string {
    return value
        .replaceAll("&", "&amp;")
        .replaceAll("<", "&lt;")
        .replaceAll(">", "&gt;")
        .replaceAll('"', "&quot;");
}
