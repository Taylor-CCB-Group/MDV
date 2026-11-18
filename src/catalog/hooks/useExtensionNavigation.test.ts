import { describe, expect, test } from "vitest";
import { parseExtensionNavigation } from "./useExtensionNavigation";

describe("parseExtensionNavigation", () => {
    test("returns valid generic navigation entries", () => {
        expect(
            parseExtensionNavigation({
                extensions: [
                    {
                        id: "example_extension",
                        label: "Example Extension",
                        url: "/example-extension/",
                    },
                ],
            }),
        ).toEqual([
            {
                id: "example_extension",
                label: "Example Extension",
                url: "/example-extension/",
            },
        ]);
    });

    test("ignores malformed entries and responses", () => {
        expect(parseExtensionNavigation({ extensions: [{ id: "incomplete" }] })).toEqual([]);
        expect(parseExtensionNavigation({})).toEqual([]);
    });
});
