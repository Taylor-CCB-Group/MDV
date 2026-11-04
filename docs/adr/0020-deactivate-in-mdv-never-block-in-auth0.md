# Removing access deactivates the User in MDV and leaves Auth0 alone

**Status:** accepted

## Context

All MDV sites share one Auth0 tenant. The deployment notes give each production site its own
password connection, such as `oceanus-cartography-username-pass`, but tell dev sites to reuse
a shared one, and a client's Microsoft connection can be enabled for more than one site.

Auth0's [block and unblock users](https://auth0.com/docs/manage-users/user-accounts/block-and-unblock-users)
page treats a block as a property of the Auth0 account that stops it "accessing applications".
Blocking a person therefore locks them out of every site on that connection. The same page
says a password reset through the Management API removes the block, and the Invitation flow
sends exactly that kind of reset (ADR 0019). A block also leaves MDV's own session in place:
Auth0's article [Does blocking a user end their session?](https://support.auth0.com/center/s/article/Does-blocking-a-user-end-their-session)
says an application's local session "will not be impacted by the block immediately".

The callback also activated any inactive User who signed in
(`auth/auth0_provider.py:400`), so a removed User could get back in by signing in.

## Decision

Removing access makes a Deactivated user in MDV and keeps their Project roles. MDV refuses
a Deactivated user at sign-in and on every later request, so a session that was open when
the Admin removed them stops working too. At sign-in they see the same no-access page as a
person with no User (ADR 0019), so the page does not show that MDV still holds their User.
MDV never blocks, unblocks or deletes the Auth0 account. The User row stays, so usage history
still names them. An Admin can reactivate them, which restores their projects.

Signing in never reactivates a User. Only an Admin's action does.

## Consequences

A Deactivated user still holds Project roles, so the share dialog and project member lists
leave them out.
