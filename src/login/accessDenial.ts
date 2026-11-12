/**
 * Why someone who signed in was not let in, as the server reports it.
 *
 * The server sends a reason code - `?denied=<reason>` on the login page - and
 * the wording lives here, with the rest of the interface. The refused account's
 * email is written into the page, never the URL: a URL ends up in browser
 * history and proxy logs.
 */

export type AccessDeniedReason = "email_unverified" | "not_set_up" | "no_account" | "setup_failed" | "unavailable";

/**
 * "retry" signs in again with the same account. "switch_account" ends the Auth0
 * session first, so the person can choose another account.
 */
export type DenialAction = {
    kind: "retry" | "switch_account";
    label: string;
};

export type AccessDenial = {
    title: string;
    detail: string;
    /** The first is the main action. */
    actions: DenialAction[];
};

const TRY_AGAIN: DenialAction = { kind: "retry", label: "Try again" };
const USE_DIFFERENT_ACCOUNT: DenialAction = { kind: "switch_account", label: "Use a different account" };

export const ACCESS_DENIED_MESSAGES: Record<AccessDeniedReason, AccessDenial> = {
    email_unverified: {
        title: "Verify your email address",
        detail: "Auth0 has not confirmed this address yet. Open the verification email Auth0 sent you and click its link, then continue.",
        actions: [{ kind: "retry", label: "I've verified my email" }, USE_DIFFERENT_ACCOUNT],
    },
    not_set_up: {
        title: "This deployment is not set up yet",
        detail: "It has no administrator yet. If you are setting it up, the server log says what is missing.",
        actions: [USE_DIFFERENT_ACCOUNT, TRY_AGAIN],
    },
    no_account: {
        title: "No access to this deployment",
        detail: "This account has not been given access here. Ask one of this deployment's administrators to add you, then try again.",
        actions: [USE_DIFFERENT_ACCOUNT, TRY_AGAIN],
    },
    setup_failed: {
        title: "Your administrator account could not be set up",
        detail: "This address is the deployment's first administrator, but setting the account up failed. The server log says why; once it is fixed, try again.",
        actions: [TRY_AGAIN, USE_DIFFERENT_ACCOUNT],
    },
    unavailable: {
        // A different account would fail the same way, so offer only a retry.
        title: "Sign-in is unavailable",
        detail: "Something went wrong while checking your account. Try again in a moment.",
        actions: [TRY_AGAIN],
    },
};

/** For a reason this page does not know - a newer server, or a hand-edited URL. */
export const GENERIC_DENIAL: AccessDenial = {
    title: "You could not be signed in",
    detail: "Try again, or use a different account. If it keeps happening, the server log says why.",
    actions: [TRY_AGAIN, USE_DIFFERENT_ACCOUNT],
};

/**
 * The account that was refused, as the server wrote it into the page - never in
 * the URL. Null when the server did not know it.
 */
export function getDeniedAccount(doc: Document = document): string | null {
    return readMeta(doc, "mdv-denied-account");
}

/**
 * The Auth0 connection the refused sign-in came through, so someone who used
 * the wrong Sign-in method can see that they did. Null when the server did not
 * know it.
 */
export function getDeniedSignInMethod(doc: Document = document): string | null {
    return readMeta(doc, "mdv-denied-sign-in-method");
}

/** Who this deployment says to ask for access, or null when it names nobody. */
export function getAccessContact(doc: Document = document): string | null {
    return readMeta(doc, "mdv-access-contact");
}

function readMeta(doc: Document, name: string): string | null {
    const content = doc.querySelector(`meta[name="${name}"]`)?.getAttribute("content");
    return content ? content : null;
}

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
