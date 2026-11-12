import { describe, expect, test } from "vitest";
import { ACCESS_DENIED_MESSAGES, GENERIC_DENIAL, getAccessDenial, getDeniedAccount } from "@/login/accessDenial";

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

describe("refusal actions", () => {
    test("every refusal offers a retry", () => {
        for (const denial of [...Object.values(ACCESS_DENIED_MESSAGES), GENERIC_DENIAL]) {
            expect(denial.actions.some((action) => action.kind === "retry")).toBe(true);
        }
    });

    test("an account with no access is pointed at switching account first", () => {
        expect(ACCESS_DENIED_MESSAGES.no_account.actions[0].kind).toBe("switch_account");
    });

    test("an outage offers only a retry, since another account would fail too", () => {
        expect(ACCESS_DENIED_MESSAGES.unavailable.actions.map((action) => action.kind)).toEqual(["retry"]);
    });
});

describe("getDeniedAccount", () => {
    function pageWithAccount(content: string | null): Document {
        const doc = document.implementation.createHTMLDocument("login");
        if (content !== null) {
            const meta = doc.createElement("meta");
            meta.name = "mdv-denied-account";
            meta.content = content;
            doc.head.appendChild(meta);
        }
        return doc;
    }

    test("reads the account the server wrote into the page", () => {
        expect(getDeniedAccount(pageWithAccount("jane@example.com"))).toBe("jane@example.com");
    });

    test("is null when the server did not know the account", () => {
        expect(getDeniedAccount(pageWithAccount(""))).toBeNull();
        expect(getDeniedAccount(pageWithAccount(null))).toBeNull();
    });
});
