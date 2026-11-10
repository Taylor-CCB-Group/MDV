/**
 * Auth0 post-login Action: tell MDV how someone signed in (ADR 0021).
 *
 * MDV accepts an Invitation only when the sign-in came through the Invitation's
 * connection and Auth0 vouches for the email. Auth0 does not return the
 * connection by default, so this Action adds it to the ID token, which also
 * puts it in the /userinfo response MDV reads. Add it to the Login flow of
 * each tenant with an MDV site.
 *
 * Auth0 runs a post-login Action for every application in the tenant, so this
 * one adds claims only for applications marked as MDV sites. Mark each MDV
 * site's application with the application metadata key "mdv" set to "true"
 * (Applications > the application > Settings > Advanced Settings >
 * Application Metadata).
 *
 * The claim names must match CLAIM_NAMESPACE in mdvtools/auth/auth0_provider.py.
 */
exports.onExecutePostLogin = async (event, api) => {
  if (event.client.metadata?.mdv !== "true") {
    return;
  }
  const namespace = "https://mdv/";
  api.idToken.setCustomClaim(`${namespace}connection`, event.connection.name);
  // "auth0" for a password connection, and "waad" or "oidc" for Microsoft
  // through Auth0's Azure AD connection or an OpenID Connect connection to
  // Entra. MDV picks the email check from it and refuses any other strategy.
  api.idToken.setCustomClaim(`${namespace}connection_strategy`, event.connection.strategy);
  if (event.connection.strategy === "waad" || event.connection.strategy === "oidc") {
    // Microsoft's claim that the tenant owns the email's domain. Only a true
    // boolean counts. Without it, MDV refuses Microsoft sign-ins that have no
    // User yet. Auth0's native Azure AD connection is reported to drop custom
    // claims, so use an OpenID Connect connection and confirm the claim
    // arrives before offering it (ADR 0021).
    //
    // Read it only from the identity this sign-in used. With linked
    // accounts, the root profile belongs to the primary identity, the first
    // in the list, and each linked identity keeps its attributes in its own
    // profileData. Nobody has confirmed where Auth0 puts the claim for a
    // primary identity, so for one read both.
    const identities = event.user.identities ?? [];
    const index = identities.findIndex((identity) => identity.connection === event.connection.name);
    const profileData = index === -1 ? {} : identities[index].profileData ?? {};
    const xmsEdov =
      profileData.xms_edov === true || (index === 0 && event.user.xms_edov === true);
    api.idToken.setCustomClaim(`${namespace}xms_edov`, xmsEdov);
  }
};
