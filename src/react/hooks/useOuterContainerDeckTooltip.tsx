import type { PickingInfo } from "@deck.gl/core";
import { useCallback, useEffect, useState, type RefObject } from "react";
import { createPortal } from "react-dom";
import { useOuterContainer } from "../screen_state";

type TooltipContent =
    | null
    | string
    | {
          text?: string;
          html?: string;
          className?: string;
          style?: Partial<CSSStyleDeclaration>;
      };

type TooltipState = {
    content: Exclude<TooltipContent, null>;
    x: number;
    y: number;
};

const TOOLTIP_OFFSET = 12;

function hasTooltipContent(content: TooltipContent | undefined): content is Exclude<TooltipContent, null> {
    if (!content) return false;
    if (typeof content === "string") return content.length > 0;
    return content.html !== undefined || content.text !== undefined;
}

function renderTooltipContent(content: Exclude<TooltipContent, null>) {
    if (typeof content === "string") return content;
    if (content.html !== undefined) {
        return <span dangerouslySetInnerHTML={{ __html: content.html }} />;
    }
    return content.text ?? null;
}

export function useOuterContainerDeckTooltip(
    getTooltipContent: (info: PickingInfo) => TooltipContent | undefined,
    anchorRef: RefObject<HTMLElement | null>,
) {
    const outerContainer = useOuterContainer();
    const [tooltip, setTooltip] = useState<TooltipState | null>(null);
    const [suppressed, setSuppressed] = useState(false);

    const clearTooltip = useCallback(() => {
        setTooltip(null);
    }, []);

    const endSuppression = useCallback(() => {
        setSuppressed(false);
    }, []);

    const suppressTooltipUntilPointerUp = useCallback(() => {
        setSuppressed(true);
        clearTooltip();
        window.removeEventListener("pointerup", endSuppression);
        window.removeEventListener("pointercancel", endSuppression);
        window.removeEventListener("blur", endSuppression);
        window.addEventListener("pointerup", endSuppression, { once: true });
        window.addEventListener("pointercancel", endSuppression, { once: true });
        window.addEventListener("blur", endSuppression, { once: true });
    }, [clearTooltip, endSuppression]);

    useEffect(
        () => () => {
            window.removeEventListener("pointerup", endSuppression);
            window.removeEventListener("pointercancel", endSuppression);
            window.removeEventListener("blur", endSuppression);
        },
        [endSuppression],
    );

    const getTooltip = useCallback(
        (info: PickingInfo) => {
            if (suppressed || !Number.isFinite(info.x) || !Number.isFinite(info.y)) {
                clearTooltip();
                return null;
            }

            const content = getTooltipContent(info);
            if (!hasTooltipContent(content)) {
                clearTooltip();
                return null;
            }

            const anchorRect = anchorRef.current?.getBoundingClientRect();
            const x = (anchorRect?.left ?? 0) + info.x + TOOLTIP_OFFSET;
            const y = (anchorRect?.top ?? 0) + info.y + TOOLTIP_OFFSET;
            setTooltip({ content, x, y });
            return null;
        },
        [anchorRef, clearTooltip, getTooltipContent, suppressed],
    );

    const tooltipPortal = tooltip
        ? createPortal(
              <div
                  className={typeof tooltip.content === "object" ? tooltip.content.className : undefined}
                  style={{
                      position: "fixed",
                      left: tooltip.x,
                      top: tooltip.y,
                      zIndex: 10000,
                      pointerEvents: "none",
                      color: "#a0a7b4",
                      backgroundColor: "#29323c",
                      padding: "10px",
                      borderRadius: 2,
                      maxWidth: 360,
                      fontSize: 12,
                      lineHeight: 1.4,
                      boxShadow: "0 4px 16px rgba(0,0,0,0.25)",
                  }}
              >
                  {renderTooltipContent(tooltip.content)}
              </div>,
              outerContainer,
          )
        : null;

    return {
        clearTooltip,
        getTooltip,
        suppressTooltipUntilPointerUp,
        tooltipPortal,
    };
}
