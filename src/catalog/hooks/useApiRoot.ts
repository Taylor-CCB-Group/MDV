import { useEffect, useState } from "react";

export function useApiRoot() {
    const [mdvApiRoot, setMdvApiRoot] = useState<string>("/");

    // Get the API root environment variable
    useEffect(() => {
        fetch("api_root")
            .then(res => res.json())
            .then(data => setMdvApiRoot(data.mdv_api_root || "/"))
            .catch(() => setMdvApiRoot("/"));
    }, []);

    return mdvApiRoot;
}