import type { PickingInfo } from "@deck.gl/core";
import { escapeHtml } from "@/utilities/Utilities";
import {
    findPickedInfo,
    getPickingInfoLayerId,
    pickingInfoMatchesLayer,
    type MDVPickingInfo,
} from "./deckPicking";

type TooltipContent =
    | null
    | string
    | {
          text?: string;
          html?: string;
          className?: string;
          style?: Partial<CSSStyleDeclaration>;
      };

type PickableInfo = PickingInfo | MDVPickingInfo | null | undefined;

const TOOLTIP_SEPARATOR = `<div style="border-top: 1px solid rgba(255,255,255,0.35); margin: 6px 0;"></div>`;

function tooltipContentToHtml(content: TooltipContent | undefined) {
    if (!content) return null;
    if (typeof content === "string") return escapeHtml(content);
    if (content.html !== undefined) return content.html;
    if (content.text !== undefined) return escapeHtml(content.text);
    return null;
}

export function combineTooltipContent(...contents: (TooltipContent | undefined)[]): TooltipContent {
    const htmlParts = contents
        .map((content) => tooltipContentToHtml(content))
        .filter((html) => html !== null);
    if (htmlParts.length === 0) return null;
    return { html: htmlParts.join(TOOLTIP_SEPARATOR) };
}

function getGateOrLabelTooltip(info: PickableInfo, gateDisplayLayerId?: string, gateLabelLayerId?: string) {
    const gateInfo = gateDisplayLayerId
        ? findPickedInfo(info, (pickedInfo) => gateDisplayLayerId === getPickingInfoLayerId(pickedInfo))
        : null;
    const gateName = gateInfo?.object?.properties?.gateName;
    if (gateName != null) {
        return {
            html: `<strong>${escapeHtml(String(gateName))}</strong><br/><small>Click on the label to edit</small>`,
        };
    }

    const labelInfo = gateLabelLayerId
        ? findPickedInfo(info, (pickedInfo) => gateLabelLayerId === getPickingInfoLayerId(pickedInfo))
        : null;
    const labelText = labelInfo?.object?.text;
    if (labelText != null) {
        return {
            html: `<strong>${escapeHtml(String(labelText))}</strong><br/><small>Click on the label to edit</small>`,
        };
    }

    return null;
}

function getScatterPickingInfo(info: PickableInfo) {
    return findPickedInfo(
        info,
        (pickedInfo) =>
            pickingInfoMatchesLayer(pickedInfo, (layerId) =>
                layerId.includes("scatter_") || layerId.includes("spatial.scatterplot"),
            ),
    );
}

export function getCombinedScatterTooltip(
    info: PickableInfo,
    options: {
        gateDisplayLayerId?: string;
        gateLabelLayerId?: string;
        getPointTooltip: (info?: PickingInfo | null) => TooltipContent | undefined;
    },
) {
    const gateOrLabelTooltip = getGateOrLabelTooltip(
        info,
        options.gateDisplayLayerId,
        options.gateLabelLayerId,
    );
    const pointTooltip = options.getPointTooltip(getScatterPickingInfo(info));
    return combineTooltipContent(gateOrLabelTooltip, pointTooltip);
}
