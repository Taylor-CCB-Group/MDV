# A first sign-in accepts an Invitation by email, connection and a verified email

**Status:** accepted

## Context

MDV cannot know a Microsoft user's Auth0 account ID before they first sign in, so the first
sign-in has to be matched to an Invitation by something the Admin typed, which is the email.
After that, MDV finds the User by their Auth0 account ID alone.

An email claim on its own is not safe to match on. Microsoft's guidance,
[Migrate away from using email claims for user identification or authorization](https://learn.microsoft.com/et-EE/azure/active-directory/develop/migrate-off-email-claim-authorization),
says the Entra email claim is mutable and unverified. The
[MSRC post on the same issue](https://msrc.microsoft.com/blog/2023/06/potential-risk-of-privilege-escalation-in-azure-ad-applications)
says apps should never use it for authorisation. On a password connection that allows
self-signup, anyone can register an account with someone else's email.

## Decision

A sign-in accepts an open Invitation only when all of these hold:

- the email matches the Invitation's email, ignoring case
- the sign-in came through the Invitation's Sign-in method, the Auth0 connection, which a
  post-login Action adds as a claim because Auth0 does not include it by default
- for the password connection, Auth0 reports the email as verified, which happens when the
  person sets their password through the emailed link
- for a Microsoft connection, the sign-in carries Microsoft's `xms_edov` claim showing that
  the email's domain is verified

Each Microsoft connection accepts sign-ins from the client's own Entra tenant only. Together
with `xms_edov`, this means only the client's tenant can produce a matching email, and only
from a domain it has proved it owns.

Both Sign-in methods use the same rule, because Microsoft Invitations will outnumber
password ones and one rule is easier to reason about.

## Consequences

Auth0's documentation does not say whether its Azure AD connection passes `xms_edov`
through. This has to be confirmed on the tenant before a Microsoft connection goes live.
If it does not arrive, Microsoft Invitations cannot be accepted until the connection or the
Action provides it.
