import { useEffect, useState } from "react";
import { apiFetch } from "@/utils/mdvRouting";

/**
 * Probes the admin session endpoint to decide whether the current user may
 * open the admin portal. The endpoint returns 200 only for admins (and for the
 * synthetic local admin when auth is disabled); it returns 401/403 otherwise,
 * so a non-ok response means the Admin entry point should stay hidden.
 */
const useAdminPortalAccess = () => {
    const [canAccessAdmin, setCanAccessAdmin] = useState(false);

    useEffect(() => {
        let cancelled = false;

        apiFetch("/admin/api/session")
            .then((res) => {
                if (!cancelled) {
                    setCanAccessAdmin(res.ok);
                }
            })
            .catch(() => {
                if (!cancelled) {
                    setCanAccessAdmin(false);
                }
            });

        return () => {
            cancelled = true;
        };
    }, []);

    return canAccessAdmin;
};

export default useAdminPortalAccess;
