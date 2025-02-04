import { useEffect, useState } from "react";

export interface User {
    name: string;
    email: string;
    association: string;
    avatarUrl: string;
}

export const mockUser: User = {
    name: "John Doe",
    email: "john.doe@example.com",
    association: "Example Corp",
    avatarUrl: "",
};

const useUser = () => {
    const [user, setUser] = useState<User | null>(null);
    const [isLoading, setIsLoading] = useState(true);
    const [error, setError] = useState<string | null>(null);

    useEffect(() => {
        const fetchUser = async () => {
            setIsLoading(true);
            try {
                // Endpoint yet to be implemented
                const response = await fetch("/user");
                if (!response.ok) {
                    throw new Error("Failed to fetch user data");
                }
                const userData: User = await response.json();

                // Check if the response is empty (you might need to adjust this condition based on your API response)
                if (Object.keys(userData).length === 0) {
                    console.warn(
                        "Received empty user data, using mock data instead.",
                    );
                    setUser(mockUser);
                } else {
                    setUser(userData);
                }
            } catch (err) {
                console.error("Error fetching user data:", err);
                console.warn("Using mock data due to error.");
                setUser(mockUser);
            } finally {
                setIsLoading(false);
            }
        };

        fetchUser();
    }, []);

    return { user, isLoading, error };
};

export default useUser;
