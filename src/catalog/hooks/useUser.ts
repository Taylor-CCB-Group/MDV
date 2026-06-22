import { useEffect, useState } from "react";
import { apiFetch } from "@/utils/mdvRouting";

export interface User {
    name: string;
    email: string;
    association: string;
    avatarUrl: string;
}

const useUser = () => {
    const [user, setUser] = useState<User | null>(null);
    const [isLoading, setIsLoading] = useState(true);
    const [error, setError] = useState<string | null>(null);

    useEffect(() => {
        const fetchUser = async () => {
            setIsLoading(true);
            setError(null);
            try {
                const response = await apiFetch("profile");
                if (!response.ok) {
                    throw new Error("Unable to fetch user data");
                }
                const userData: User = await response.json();

                // Check if the response is empty
                if (Object.keys(userData).length === 0) {
                    throw new Error("Received empty user data");
                }
                setUser(userData);
            } catch (err) {
                console.error("Error fetching user data:", err);
                setUser(null);
                setError("Unable to fetch user data");
            } finally {
                setIsLoading(false);
            }
        };

        fetchUser();
    }, []);

    return { user, isLoading, error };
};

export default useUser;
