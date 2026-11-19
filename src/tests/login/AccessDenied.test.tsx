import { render, screen } from "@testing-library/react";
import { describe, expect, test } from "vitest";
import AccessDenied from "@/login/AccessDenied";
import { ACCESS_DENIED_MESSAGES } from "@/login/accessDenial";

// ADR 0019: a refused person sees the email and Sign-in method they used, and
// the site contact when one is configured.

describe("AccessDenied", () => {
    test("names the account, how it signed in and who to ask", () => {
        render(
            <AccessDenied
                denial={ACCESS_DENIED_MESSAGES.no_account}
                account="ann@acme.example"
                signInMethod="acme-entra"
                contact="MDV support, mdv-support@example.org"
            />,
        );

        expect(screen.getByText("ann@acme.example")).toBeDefined();
        expect(screen.getByText("acme-entra")).toBeDefined();
        expect(screen.getByText(/contact MDV support, mdv-support@example.org/)).toBeDefined();
    });

    test("leaves out what the server did not know", () => {
        render(
            <AccessDenied
                denial={ACCESS_DENIED_MESSAGES.no_account}
                account="ann@acme.example"
                signInMethod={null}
                contact={null}
            />,
        );

        expect(screen.getByText("ann@acme.example")).toBeDefined();
        expect(screen.queryByText(/through/)).toBeNull();
        expect(screen.queryByText(/contact/)).toBeNull();
    });
});
