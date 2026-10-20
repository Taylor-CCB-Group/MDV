/**
 * Why someone who signed in was not let in, as the server reports it.
 *
 * The server only sends a reason code - `?denied=<reason>` on the login page -
 * and the wording lives here, with the rest of the interface. It sends no
 * address: an email in a URL ends up in browser history and proxy logs.
 */

export type AccessDeniedReason = "email_unverified" | "not_set_up" | "no_account" | "setup_failed" | "unavailable";

export type AccessDenial = {
    title: string;
    detail: string;
};

export const ACCESS_DENIED_MESSAGES: Record<AccessDeniedReason, AccessDenial> = {
    email_unverified: {
        title: "Verify your email address",
        detail: "Auth0 has not confirmed this address yet. Open the verification email Auth0 sent you, click its link, then try again.",
    },
    not_set_up: {
        title: "This deployment is not set up yet",
        detail: "It has no administrator yet. If you are setting it up, the server log says what is missing.",
    },
    no_account: {
        title: "No access to this deployment",
        detail: "This account has not been given access here. Ask one of this deployment's administrators to add you, or sign out and use another account.",
    },
    setup_failed: {
        title: "Your administrator account could not be set up",
        detail: "This address is the deployment's first administrator, but setting the account up failed. The server log says why; once it is fixed, try again.",
    },
    unavailable: {
        title: "Sign-in is unavailable",
        detail: "Something went wrong while checking your account. Try again in a moment.",
    },
};

/** For a reason this page does not know - a newer server, or a hand-edited URL. */
export const GENERIC_DENIAL: AccessDenial = {
    title: "You could not be signed in",
    detail: "Try again, or sign out and use another account. If it keeps happening, the server log says why.",
};

function isAccessDeniedReason(value: string): value is AccessDeniedReason {
    return Object.hasOwn(ACCESS_DENIED_MESSAGES, value);
}

/**
 * The denial to show, or null when the login page should show sign-in.
 *
 * An unknown reason still gets a page rather than sign-in: offering a sign-in
 * button to someone the server just turned away is how the old loop started.
 */
export function getAccessDenial(search: string = window.location.search): AccessDenial | null {
    const reason = new URLSearchParams(search).get("denied");
    if (reason === null) return null;
    return isAccessDeniedReason(reason) ? ACCESS_DENIED_MESSAGES[reason] : GENERIC_DENIAL;
}
