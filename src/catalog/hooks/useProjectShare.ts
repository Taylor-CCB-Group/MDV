import { useEffect, useState } from "react";

export type UserPermission = "View" | "Edit" | "Owner";

export type SharedUser = {
    id: number;
    email: string;
    permission: UserPermission;
};

export type RegisteredUser = {
    id: number;
    email: string;
};

const useProjectShare = (projectId: string) => {
    const [email, setEmail] = useState("");
    const [sharedUsers, setSharedUsers] = useState<SharedUser[]>([]);
    const [newUser, setNewUser] = useState<RegisteredUser | null>();
    const [userList, setUserList] = useState<RegisteredUser[]>([]);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState("");
    const [errorMsg, setErrorMsg] = useState("");

    useEffect(() => {
        if (projectId) {
            getAllUsers();
        }
    }, [projectId]);

    const getAllUsers = async () => {
        setError("");
        setErrorMsg("");
        setIsLoading(true);
        try {
            const res = await fetch(`projects/${projectId}/share`, {
                headers: {
                    Accept: "application/json",
                },
            });

            const data = await res.json();

            if (data?.error) {
                setError(data?.error);
                return;
            }

            if (data?.status === "error") {
                setErrorMsg(data?.message);
                return;
            }

            if (res.ok) {
                if (data?.all_users) setUserList(data.all_users);
                if (data?.shared_users) setSharedUsers(data.shared_users);
                return;
            } else {
                if (res.status === 403) {
                    setError(data?.error || "Forbidden: You are not allowed to perform this action.");
                } else {
                    setError(data?.error || "Internal Server Error");
                }
            }

            setError(data?.message || "An error occurred. Please try again.");
        } catch (error) {
            setError("Error fetching user data. Please try again.");
        } finally {
            setIsLoading(false);
        }
    };

    const addUser = async (userId: number, permission: UserPermission) => {
        setError("");
        setErrorMsg("");
        setIsLoading(true);
        try {
            const res = await fetch(`projects/${projectId}/share`, {
                method: "POST",
                body: JSON.stringify({ user_id: userId, permission: permission.toLowerCase() }),
                headers: {
                    "Content-Type": "application/json",
                    Accept: "application/json",
                },
            });

            const data = await res.json();

            if (res.ok) {
                await getAllUsers();
            } else {
                if (res.status === 403) {
                    setError(data?.error || "Forbidden: You are not allowed to perform this action.");
                } else {
                    setError(data?.error || "Internal Server Error");
                }
            }

            setErrorMsg(data?.message || "An error occurred. Please try again.");
        } catch (error) {
            setErrorMsg("Error adding new user. Please try again.");
        } finally {
            setIsLoading(false);
        }
    };

    const changeUserPermission = async (userId: number, permission: UserPermission) => {
        setError("");
        setErrorMsg("");
        setIsLoading(true);
        try {
            const res = await fetch(`projects/${projectId}/share/${userId}/edit`, {
                method: "POST",
                body: JSON.stringify({ permission }),
                headers: {
                    "Content-Type": "application/json",
                    Accept: "application/json",
                },
            });

            const data = await res.json();

            if (res.ok) {
                await getAllUsers();
            } else {
                if (res.status === 403) {
                    setError(data?.error || "Forbidden: You are not allowed to perform this action.");
                } else {
                    setError(data?.error || "Internal Server Error");
                }
            }

            setErrorMsg(data?.message || "An error occurred. Please try again.");
        } catch (error) {
            setErrorMsg("Error fetching user data. Please try again.");
        } finally {
            setIsLoading(false);
        }
    };

    const deleteSharedUser = async (userId: number) => {
        setError("");
        setErrorMsg("");
        setIsLoading(true);
        try {
            const res = await fetch(`projects/${projectId}/share/${userId}/delete`, {
                method: "POST",
                headers: {
                    Accept: "application/json",
                },
            });

            const data = await res.json();

            if (res.ok) {
                await getAllUsers();
            } else {
                if (res.status === 403) {
                    setError(data?.error || "Forbidden: You are not allowed to perform this action.");
                } else {
                    setError(data?.error || "Internal Server Error");
                }
            }

            setErrorMsg(data?.message || "An error occurred. Please try again.");
        } catch (error) {
            setErrorMsg("Error deleting user. Please try again.");
        } finally {
            setIsLoading(false);
        }
    };

    return {
        email,
        setEmail,
        sharedUsers,
        setSharedUsers,
        addUser,
        getAllUsers,
        isLoading,
        error,
        errorMsg,
        userList,
        newUser,
        setNewUser,
        changeUserPermission,
        deleteSharedUser,
    };
};

export default useProjectShare;
