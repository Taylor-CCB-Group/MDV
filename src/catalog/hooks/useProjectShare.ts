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


            if (res.ok) {
                const data = await res.json();
                if (data?.all_users) setUserList(data.all_users);
                if (data?.shared_users) setSharedUsers(data.shared_users);
                return;
            } else {
                const errorData = await res.json().catch(() => ({
                    error: "An error occurred. Please try again."
                }))
                if (res.status === 403) {
                    setError(errorData?.error || "Forbidden: You are not allowed to perform this action.");
                } else if (res.status === 500) {
                    setError(errorData?.error || "Internal Server Error");
                } else {
                    setError(errorData?.error || "An error occurred. Please try again.");
                }
            }

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


            if (res.ok) {
                const data = await res.json();
                await getAllUsers();
                return;
            } else {
                const errorData = await res.json().catch(() => ({
                    error: "An error occurred. Please try again."
                }))
                if (res.status === 403) {
                    setError(errorData?.error || "Forbidden: You are not allowed to perform this action.");
                } else if (res.status === 500) {
                    setError(errorData?.error || "Internal Server Error");
                } else {
                    setErrorMsg(errorData?.message || "An error occurred. Please try again.");
                }
            }

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


            if (res.ok) {
                const data = await res.json();
                await getAllUsers();
                return;
            } else {
                const errorData = await res.json().catch(() => ({
                    error: "An error occurred. Please try again."
                }))
                if (res.status === 403) {
                    setError(errorData?.error || "Forbidden: You are not allowed to perform this action.");
                } else if (res.status === 500) {
                    setError(errorData?.error || "Internal Server Error");
                } else {
                    setErrorMsg(errorData?.message || "An error occurred. Please try again.");
                }
            }

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


            if (res.ok) {
                const data = await res.json();
                await getAllUsers();
                return;
            } else {
                const errorData = await res.json().catch(() => ({
                    error: "An error occurred. Please try again."
                }))
                if (res.status === 403) {
                    setError(errorData?.error || "Forbidden: You are not allowed to perform this action.");
                } else {
                    setError(errorData?.error || "Internal Server Error");
                }
                setErrorMsg(errorData?.message || "An error occurred. Please try again.");
            }

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
