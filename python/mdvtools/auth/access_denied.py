"""Why a signed-in person was turned away.

Someone can sign in to Auth0 and still not be allowed into this deployment - no
account here, an unverified email, a failed first-administrator setup. Sending
them to the login page as if they were simply signed out does not work: Auth0
still has their session, signs them straight back in, and they land on the same
refusal again, in a loop.

So the provider records a reason during the request, and the auth hook sends
the person to the login page with that reason in the URL. The page shows what
went wrong instead of a sign-in button. The wording lives in the frontend, in
src/login/accessDenial.ts; only the reason code travels.
"""

from typing import Literal, Optional

from flask import g, session

AccessDeniedReason = Literal["email_unverified", "not_set_up", "no_account", "setup_failed", "unavailable"]

# Read by templates/login.html, which names the refused account on the page.
DENIED_ACCOUNT_SESSION_KEY = "access_denied_email"
# The Auth0 connection the refused sign-in came through, so someone who used the
# wrong Sign-in method can see that they did (ADR 0019).
DENIED_SIGN_IN_METHOD_SESSION_KEY = "access_denied_sign_in_method"


def deny_access(
    reason: AccessDeniedReason, email: Optional[str] = None, sign_in_method: Optional[str] = None
) -> None:
    """Record why this request's signed-in user was refused.

    The reason is stored on flask.g, which lasts for one request only, so the
    next sign-in attempt is judged afresh - that is what makes "verify your
    email, then try again" work. The email and Sign-in method, when known, go
    in the session rather than the URL, so the login page can say which account
    was refused and how it signed in.
    """
    g.access_denied_reason = reason
    for key, value in ((DENIED_ACCOUNT_SESSION_KEY, email), (DENIED_SIGN_IN_METHOD_SESSION_KEY, sign_in_method)):
        if value:
            session[key] = value
        else:
            session.pop(key, None)


def access_denied_reason() -> Optional[AccessDeniedReason]:
    """The reason recorded by deny_access during this request, if any."""
    return g.get("access_denied_reason")


def login_url_with_reason(login_url: str, reason: AccessDeniedReason) -> str:
    """The login page URL with ?denied=<reason> on the end.

    The reason is one of the fixed codes above, so it needs no escaping. The
    email address is deliberately not included: a URL ends up in browser
    history and proxy logs.
    """
    separator = "&" if "?" in login_url else "?"
    return f"{login_url}{separator}denied={reason}"
