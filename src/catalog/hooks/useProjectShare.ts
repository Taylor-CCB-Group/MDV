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

const users: SharedUser[] = [
    { id: 1, email: "joe@test.com", permission: "Edit" },
    { id: 2, email: "pete@test.com", permission: "View" },
    { id: 3, email: "matt@test.com", permission: "View" },
];

const mockUserList: RegisteredUser[] = [
    { id: 1, email: "joe@test.com" },
    { id: 2, email: "pete@test.com" },
    { id: 3, email: "matt@test.com" },
];

const useProjectShare = (projectId: string) => {
    const [email, setEmail] = useState("");
    const [sharedUsers, setSharedUsers] = useState<SharedUser[]>([]);
    const [newUser, setNewUser] = useState<RegisteredUser>();
    const [userList, setUserList] = useState<RegisteredUser[]>([]);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState<Error>();

    useEffect(() => {
        if (projectId) {
            // setUserList(mockUserList);
            // setSharedUsers(users);
            // todo: uncomment later
            getAllUsers();
        }
    }, [projectId]);

    const getAllUsers = async () => {
        setIsLoading(true);
        try {
            // todo: Change the api endpoint
            const res = await fetch(`projects/${projectId}/share`, {
                headers: {
                    Accept: "application/json",
                },
            });
            if (res.ok) {
                const data = await res.json();
                // todo: update user list with this data
                if (data?.all_users)    setUserList(data.all_users);
                if (data?.shared_users) setSharedUsers(data.shared_users);
                console.log("getAllUsers", data);
            } else {
                throw new Error("An unknown error occurred.")
            }


            if (res.status === 500) {
                throw new Error("Internal Server Error.");
            }
        } catch (error) {
            const err = error instanceof Error ? error : new Error("Error fetching users data.");
            setError(err);
        } finally {
            setIsLoading(false);
        }
    };

    const addUser = async (userId: number, permission: UserPermission) => {
        setIsLoading(true);
        try {
            // todo: Change the api endpoint
            const res = await fetch(`projects/${projectId}/share`, {
                method: "POST",
                body: JSON.stringify({ userId, permission }),
            });

            console.log("addUser", res);

            if (!res.ok) {
                throw new Error("An unknown error occurred.");
            }

            if (res.status === 500) {
                throw new Error("Internal Server Error.");
            }

            // todo: update shared users
            // setSharedUsers();
            await getAllUsers();

        } catch (error) {
            const err = error instanceof Error ? error : new Error("Error adding new user.");
            setError(err);
        } finally {
            setIsLoading(false);
        }
    };

    const changeUserPermission = async (userId: number, permission: UserPermission) => {
        setIsLoading(true);
        try {
            // todo: Change the api endpoint
            const res = await fetch(`/projects/${projectId}/share/${userId}/edit`, {
                method: "POST",
                body: JSON.stringify({ userId, permission }),
            });

            console.log("changeUserPermission", res);

            if (!res.ok) {
                throw new Error("An unknown error occurred.");
            }

            if (res.status === 500) {
                throw new Error("Internal Server Error.");
            }

            // todo: update shared users
            // setSharedUsers();
            await getAllUsers();
        } catch (error) {
            const err = error instanceof Error ? error : new Error("Error changing user permission.");
            setError(err);
        } finally {
            setIsLoading(false);
        }
    };

    const deleteSharedUser = async (userId: number) => {
        setIsLoading(true);
        try {
            // todo: Change the api endpoint
            const res = await fetch(`/projects/${projectId}/share/${userId}/edit`, {
                method: "POST",
                // body: JSON.stringify({ user }),
            });

            console.log("deleteSharedUser", res);

            if (!res.ok) {
                throw new Error("An unknown error occurred.");
            }

            if (res.status === 500) {
                throw new Error("Internal Server Error.");
            }

            // todo: update shared users
            // setSharedUsers();
            await getAllUsers();

        } catch (error) {
            const err = error instanceof Error ? error : new Error("Error deleting user.");
            setError(err);
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
        userList,
        newUser,
        setNewUser,
        changeUserPermission,
        deleteSharedUser,
    };
};

export default useProjectShare;
