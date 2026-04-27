import type { PickingInfo } from "@deck.gl/core";

type DeckPicker = {
    pickMultipleObjects: (options: {
        x: number;
        y: number;
        radius?: number;
        depth?: number;
        layerIds?: string[];
        unproject3D?: boolean;
    }) => PickingInfo[];
};

export type MDVPickingInfo = PickingInfo & {
    mdvPickedInfos: PickingInfo[];
};

function isDeckPicker(deck: unknown): deck is DeckPicker {
    return (
        typeof deck === "object" &&
        deck !== null &&
        "pickMultipleObjects" in deck &&
        typeof deck.pickMultipleObjects === "function"
    );
}

function hasMDVPickedInfos(info: PickingInfo | MDVPickingInfo | null | undefined): info is MDVPickingInfo {
    return !!info && "mdvPickedInfos" in info && Array.isArray(info.mdvPickedInfos);
}

function getUniquePickedInfos(primary: PickingInfo, pickedInfos: PickingInfo[]) {
    const uniqueInfos: PickingInfo[] = [];
    const keys = new Set<string>();

    for (const info of [primary, ...pickedInfos]) {
        const layerId = getPickingInfoLayerId(info);
        const sourceLayerId = getPickingInfoSourceLayerId(info);
        const key = `${layerId}|${sourceLayerId}|${info.index}`;
        if (keys.has(key)) continue;
        keys.add(key);
        uniqueInfos.push(info);
    }

    return uniqueInfos;
}

export function getPickingInfoWithAlternates(
    info: PickingInfo,
    deck: unknown,
    options?: { depth?: number; radius?: number },
): MDVPickingInfo;
export function getPickingInfoWithAlternates(
    info: PickingInfo | null | undefined,
    deck: unknown,
    options?: { depth?: number; radius?: number },
): MDVPickingInfo | null;
export function getPickingInfoWithAlternates(
    info: PickingInfo | null | undefined,
    deck: unknown,
    options: { depth?: number; radius?: number } = {},
): MDVPickingInfo | null {
    if (!info) return null;

    const { x, y } = info;
    if (!isDeckPicker(deck) || !Number.isFinite(x) || !Number.isFinite(y)) {
        return { ...info, mdvPickedInfos: [info] };
    }

    const pickedInfos = deck.pickMultipleObjects({
        x,
        y,
        radius: options.radius ?? 2,
        depth: options.depth ?? 8,
    });

    return {
        ...info,
        mdvPickedInfos: getUniquePickedInfos(info, pickedInfos),
    };
}

export function getPickedInfos(info: PickingInfo | MDVPickingInfo | null | undefined) {
    if (!info) return [];
    return hasMDVPickedInfos(info) ? info.mdvPickedInfos : [info];
}

export function findPickedInfo(
    info: PickingInfo | MDVPickingInfo | null | undefined,
    predicate: (info: PickingInfo) => boolean,
) {
    return getPickedInfos(info).find(predicate) ?? null;
}

export function getPickingInfoLayerId(info: PickingInfo | null | undefined) {
    return info?.layer?.id ?? "";
}

export function getPickingInfoSourceLayerId(info: PickingInfo | null | undefined) {
    return info?.sourceLayer?.id ?? "";
}

export function pickingInfoMatchesLayer(
    info: PickingInfo | null | undefined,
    predicate: (layerId: string) => boolean,
) {
    const layerId = getPickingInfoLayerId(info);
    const sourceLayerId = getPickingInfoSourceLayerId(info);
    return predicate(layerId) || predicate(sourceLayerId);
}
