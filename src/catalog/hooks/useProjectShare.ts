import axios from "axios";
import { useEffect, useState } from "react";

export type SharedUser = {
    name: string;
    permission: string;
};

const users = [
    { name: "Joe", permission: "edit" },
    { name: "Pete", permission: "view" },
    { name: "Matt", permission: "view" },
];

const mockUserList = ["Joe", "Pete", "Matt"];

const useProjectShare = (projectId: string) => {
    const [username, setUsername] = useState("");
    const [sharedUsers, setSharedUsers] = useState<SharedUser[]>([]);
    const [userList, setUserList] = useState<string[]>([]);
    const [isLoading, setIsLoading] = useState(false);
    const [errorMessage, setErrorMessage] = useState("");

    useEffect(() => {
        if (projectId) {
            setUserList(mockUserList);
            setSharedUsers(users);
            // todo: uncomment later
            // getAllUsers();
        }
    }, [projectId]);

    const getAllUsers = async () => {
        setIsLoading(true);
        try {
            // todo: Change the api endpoint
            const res = await fetch(`users/${projectId}`, {
                headers: {
                    Accept: "application/json",
                },
            });
            if (!res.ok) {
                throw new Error(`HTTP error! status: ${res.status}`);
            }
            const data = await res.json();
            console.log("getAllUsers", res);
        } catch (error) {
            console.log("error", error);
        } finally {
            setIsLoading(false);
        }
    };

    const addUser = async (username: string) => {
        setIsLoading(true);
        try {
            // todo: Change the api endpoint
            const res = await fetch("add_user", {
                method: "POST",
                body: JSON.stringify({ username }),
            });

            console.log("addUser", res);
        } catch (error) {
        } finally {
            setIsLoading(false);
        }
    };

    const updateSharedUsers = async (updatedUsers: any) => {
        setIsLoading(true);
        try {
            // todo: Change the api endpoint
            const res = await fetch("update_users", {
                method: "POST",
                body: JSON.stringify({ users: updatedUsers }),
            });

            console.log("updateSharedUsers", res);
        } catch (error) {
        } finally {
            setIsLoading(false);
        }
    };

    return {
        username,
        setUsername,
        sharedUsers,
        setSharedUsers,
        addUser,
        updateSharedUsers,
        getAllUsers,
        isLoading,
        errorMessage,
        userList
    };
};

export default useProjectShare;
