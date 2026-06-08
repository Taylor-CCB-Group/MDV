import type { Device } from "@luma.gl/core";

/** Resolve luma device from a deck.gl Deck instance (internal layout varies by version). */
export function resolveDeckDevice(deck: unknown): Device | null {
    if (!deck || typeof deck !== "object") return null;
    const typed = deck as Record<string, unknown>;
    if (typed.device) return typed.device as Device;
    const animationLoop = typed.animationLoop as Record<string, unknown> | undefined;
    if (animationLoop?.device) return animationLoop.device as Device;
    const layerManager = typed.layerManager as Record<string, unknown> | undefined;
    const context = layerManager?.context as Record<string, unknown> | undefined;
    if (context?.device) return context.device as Device;
    return null;
}
