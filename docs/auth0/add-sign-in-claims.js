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
    // Microsoft's claim that the tenant owns the email's domain. Nobody has
    // confirmed yet whether Auth0 puts it on the identity's profileData or on
    // event.user, so read both. Only a true boolean counts. Without it, MDV
    // refuses Microsoft sign-ins that have no User yet.
    const profileData = event.user.identities?.[0]?.profileData ?? {};
    const xmsEdov = profileData.xms_edov === true || event.user.xms_edov === true;
    api.idToken.setCustomClaim(`${namespace}xms_edov`, xmsEdov);
  }
};
