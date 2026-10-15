import { ACCESS_DENIED_MESSAGES, GENERIC_DENIAL, getAccessDenial } from "@/login/accessDenial";
import { describe, expect, test } from "vitest";

describe("getAccessDenial", () => {
    test("shows sign-in when the server sent no reason", () => {
        expect(getAccessDenial("")).toBeNull();
        expect(getAccessDenial("?callback=1")).toBeNull();
    });

    test("shows the message for each reason the server sends", () => {
        // The codes the Python side passes to deny_access(); a reason missing
        // here would fall back to the generic message instead of its own.
        for (const reason of ["email_unverified", "not_set_up", "no_account", "setup_failed", "unavailable"] as const) {
            expect(getAccessDenial(`?denied=${reason}`)).toBe(ACCESS_DENIED_MESSAGES[reason]);
        }
    });

    test("finds the reason among other parameters", () => {
        expect(getAccessDenial("?callback=1&denied=no_account")).toBe(ACCESS_DENIED_MESSAGES.no_account);
    });

    test("still refuses, generically, for a reason it does not know", () => {
        // Falling back to sign-in here would restart the loop.
        expect(getAccessDenial("?denied=something_new")).toBe(GENERIC_DENIAL);
        expect(getAccessDenial("?denied=")).toBe(GENERIC_DENIAL);
    });

    test("does not treat inherited object keys as reasons", () => {
        expect(getAccessDenial("?denied=toString")).toBe(GENERIC_DENIAL);
        expect(getAccessDenial("?denied=__proto__")).toBe(GENERIC_DENIAL);
    });
});
