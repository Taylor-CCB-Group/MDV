import { defineChoiceSetting } from "@/lib/browserSettings";
import { autorun } from "mobx";
import { afterEach, describe, expect, test } from "vitest";

const KEY = "mdv.settings.test";

function define(legacy?: () => "a" | "b" | undefined) {
    return defineChoiceSetting<"a" | "b">({
        key: KEY,
        label: "Test",
        choices: [
            { value: "a", label: "A" },
            { value: "b", label: "B" },
        ],
        defaultValue: "a",
        legacy,
    });
}

afterEach(() => localStorage.clear());

describe("defineChoiceSetting", () => {
    test("uses the default when nothing is stored", () => {
        expect(define().get()).toBe("a");
    });

    test("persists a set value and reads it back", () => {
        define().set("b");
        expect(localStorage.getItem(KEY)).toBe("b");
        expect(define().get()).toBe("b");
    });

    test("ignores a stored value that is not a choice", () => {
        localStorage.setItem(KEY, "nonsense");
        expect(define().get()).toBe("a");
    });

    test("reads the legacy value only when nothing current is stored", () => {
        expect(define(() => "b").get()).toBe("b");
        localStorage.setItem(KEY, "a");
        expect(define(() => "b").get()).toBe("a");
    });

    test("is observable", () => {
        const setting = define();
        const seen: string[] = [];
        const dispose = autorun(() => seen.push(setting.get()));
        setting.set("b");
        dispose();
        expect(seen).toEqual(["a", "b"]);
    });
});
