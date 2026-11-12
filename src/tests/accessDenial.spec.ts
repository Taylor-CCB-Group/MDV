import { describe, expect, test } from "vitest";
import {
    ACCESS_DENIED_MESSAGES,
    GENERIC_DENIAL,
    getAccessContact,
    getAccessDenial,
    getDeniedAccount,
    getDeniedSignInMethod,
} from "@/login/accessDenial";

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

describe("what else the server wrote into the page", () => {
    // ADR 0019: someone refused sees how they signed in and who to ask.
    function pageWith(meta: Record<string, string>): Document {
        const doc = document.implementation.createHTMLDocument("login");
        for (const [name, content] of Object.entries(meta)) {
            const tag = doc.createElement("meta");
            tag.name = name;
            tag.content = content;
            doc.head.appendChild(tag);
        }
        return doc;
    }

    test("reads the Sign-in method and the site contact", () => {
        const doc = pageWith({
            "mdv-denied-sign-in-method": "acme-entra",
            "mdv-access-contact": "MDV support, mdv-support@example.org",
        });
        expect(getDeniedSignInMethod(doc)).toBe("acme-entra");
        expect(getAccessContact(doc)).toBe("MDV support, mdv-support@example.org");
    });

    test("is null for each when the server left it empty or out", () => {
        const doc = pageWith({ "mdv-denied-sign-in-method": "" });
        expect(getDeniedSignInMethod(doc)).toBeNull();
        expect(getAccessContact(doc)).toBeNull();
    });
});
