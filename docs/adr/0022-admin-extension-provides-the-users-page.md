# An admin extension provides the users page, and MDV core keeps sign-in and access

**Status:** accepted

## Context

ADR 0019 gives Admins a users page for Invitations, Deactivated users and Admin rights. MDV
loads installed packages that register under the `mdv.extensions` entry point and are listed in
the deployment's `extensions` setting (`dbutils/external_extensions.py`). An admin extension
built on that interface already provides a users page, password invitations through Auth0 and
per-project access controls. It runs inside MDV's Flask app and uses MDV's models, and it owns
no tables.

Deployments with `ENABLE_AUTH` off, such as local installs, have no Users to manage.

## Decision

The admin extension provides the users page and everything an Admin does on it. That covers
sending, resending and revoking Invitations, deactivating and reactivating Users, and giving or
removing Admin rights. Every deployment with `ENABLE_AUTH` on installs and enables it, and
deployments with auth off leave it out.

MDV core keeps everything that runs when a person signs in or makes a request. That covers the
Invitation table with its expiry and email helpers, matching an Invitation at the callback,
refusing a Deactivated user on every request, and the no-access page.

The admin extension makes its changes by writing MDV's models directly, as it already does for
Users and Project roles. Core offers no separate service for them.

## Consequences

MDV's published image does not include the admin extension, so auth-enabled deployments run an
image that adds it. An auth-enabled deployment that does not enable the extension still starts,
but once the sync is removed nobody on it can invite anyone.

## Considered options

- **Build the users page into MDV core.** Rejected so that admin features stay an add-on that
  each deployment chooses to install.
- **Give core a service that the extension calls for every change.** Rejected because most of
  its functions would have no caller in core, and every auth-enabled deployment installs the
  extension alongside core.
