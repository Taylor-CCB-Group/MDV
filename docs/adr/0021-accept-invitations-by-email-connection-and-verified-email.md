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
- the sign-in came through the Invitation's Sign-in method, the Auth0 connection
- for the password connection, Auth0 reports the email as verified
- for a Microsoft connection, the sign-in carries Microsoft's `xms_edov` claim showing that
  the email's domain is verified

Auth0 does not include the connection in what it returns by default, so a post-login Action
adds the connection's name and strategy as claims. The strategy says which of the last 2
checks applies. A Microsoft connection is either Auth0's native Azure AD connection
(strategy `waad`) or an Enterprise OpenID Connect connection to Entra (strategy `oidc`).
MDV refuses a sign-in through a connection with any other strategy.

MDV-Admin creates a password account marked as verified, with a random password that it
discards. Only someone who can read that inbox can then set a password, through the reset
email Auth0 sends there. MDV-Admin does not change an account that existed before the
Invitation, so if that account's email is unverified, its owner has to verify it through
Auth0 before the Invitation can be accepted.

Each Microsoft connection accepts sign-ins from the client's own Entra tenant only. Together
with `xms_edov`, this means only the client's tenant can produce a matching email, and only
from a domain it has proved it owns.

Both Sign-in methods use the same rule, because Microsoft Invitations will outnumber
password ones and one rule is easier to reason about. A Deactivated user who signs in with
the Auth0 account tied to their User accepts an Invitation that names them by the same
rule, and accepting it reactivates them.

## Consequences

Auth0's documentation says its native Azure AD connection ignores Entra's optional claims,
while an Enterprise OpenID Connect connection passes them through, so `xms_edov` probably
needs an OpenID Connect connection.

Someone with rights over the client's Entra app registration has to request `xms_edov`.
[Microsoft's optional claims reference](https://learn.microsoft.com/en-us/entra/identity-platform/optional-claims-reference)
says Entra returns it only when the token also carries the `email` claim. The Azure portal
does not offer `xms_edov` in its list, so
[Clerk's guide to Microsoft sign-in](https://clerk.com/docs/authentication/social-connections/microsoft)
adds it through the manifest:

1. Under Token configuration, add the optional claims `email` and `xms_pdl` for the ID
   token, and turn on the Microsoft Graph email permission when the portal asks. Repeat
   this for the access token.
2. In the Manifest, set `acceptMappedClaims` to `true`.
3. In the manifest's `optionalClaims`, rename each `xms_pdl` entry to `xms_edov`, then save.

The portal may then warn that `xms_edov` is not recognised. Clerk's guide says the warning
can be ignored. Microsoft's reference lists the claim but does not describe these steps.

The post-login Action reads the claim only from the identity the sign-in used, so a linked
account cannot pass on another connection's verification. Auth0's
[account linking documentation](https://auth0.com/docs/manage-users/user-accounts/user-account-linking)
keeps the primary identity's attributes on the root profile and puts each linked identity's
attributes in its own `profileData`. Nobody has confirmed where the claim lands for a
primary identity, so for one the Action reads both. The connection type, the claim and where
it lands have to be confirmed on the tenant before a Microsoft connection goes live, for
example by logging `event.user` from a test Action on that connection. If the
claim does not arrive, Microsoft Invitations cannot be accepted until the connection or the
Action provides it.
