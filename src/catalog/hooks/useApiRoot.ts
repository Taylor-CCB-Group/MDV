import { useEffect, useState } from "react";

export function useApiRoot() {
    const [mdvApiRoot, setMdvApiRoot] = useState<string>("/");

    // Get the API root environment variable
    useEffect(() => {
        //! This is a relative route, we need to change it at some point to avoid issues when in a nested route
        fetch("api_root")
            .then(res => res.json())
            .then(data => setMdvApiRoot(data.mdv_api_root || "/"))
            .catch(() => setMdvApiRoot("/"));
    }, []);

    return mdvApiRoot;
}