import { useEffect, useState } from "react";
import { getDashboardApiRoot } from "@/utils/mdvRouting";

export function useApiRoot() {
    const [mdvApiRoot, setMdvApiRoot] = useState<string>(getDashboardApiRoot());

    useEffect(() => {
        setMdvApiRoot(getDashboardApiRoot());
    }, []);

    return mdvApiRoot;
}
