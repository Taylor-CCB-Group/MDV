import { useEffect, useState } from "react";
import { apiFetch, getDashboardApiRoot } from "@/utils/mdvRouting";

export function useApiRoot() {
    const [mdvApiRoot, setMdvApiRoot] = useState<string>(getDashboardApiRoot());

    useEffect(() => {
        apiFetch("api_root")
            .then((res) => res.json())
            .then((data) => setMdvApiRoot(data.mdv_api_root || getDashboardApiRoot()))
            .catch(() => setMdvApiRoot(getDashboardApiRoot()));
    }, []);

    return mdvApiRoot;
}
