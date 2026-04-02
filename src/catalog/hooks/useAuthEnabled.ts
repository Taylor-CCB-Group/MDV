import { useEffect, useState } from "react";
import { apiFetch } from "@/utils/mdvRouting";

const useAuthEnabled = () => {
    const [authEnabled, setAuthEnabled] = useState(false);
    useEffect(() => {
        apiFetch("enable_auth")
            .then(res => res.json())
            .then(data => setAuthEnabled(data.enable_auth || false))
            .catch(() => setAuthEnabled(false));
    }, []);

    return authEnabled;
};

export default useAuthEnabled;
