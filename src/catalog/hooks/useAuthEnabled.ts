import { useEffect, useState } from "react";

const useAuthEnabled = () => {
    const [authEnabled, setAuthEnabled] = useState(false);
    useEffect(() => {
        fetch("enable_auth")
            .then(res => res.json())
            .then(data => setAuthEnabled(data.enable_auth || false))
            .catch(() => setAuthEnabled(false));
    }, []);

    return authEnabled;
};

export default useAuthEnabled;