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

from flask import g

AccessDeniedReason = Literal["email_unverified", "not_set_up", "no_account", "setup_failed", "unavailable"]


def deny_access(reason: AccessDeniedReason) -> None:
    """Record why this request's signed-in user was refused.

    Stored on flask.g, which lasts for one request only, so the next sign-in
    attempt is judged afresh - that is what makes "verify your email, then try
    again" work.
    """
    g.access_denied_reason = reason


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
