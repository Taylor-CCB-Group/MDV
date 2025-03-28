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

const useProjectShare = (projectId: string) => {
    const [username, setUsername] = useState("");
    const [sharedUsers, setSharedUsers] = useState<SharedUser[]>([]);

    useEffect(() => {
        if (projectId) {
            setSharedUsers(users);
            // todo: uncomment later
            // getAllUsers();
        }
    }, [projectId]);

    const getAllUsers = async () => {
        // todo: Change the api endpoint
        const res = await fetch(`users/${projectId}`, {
            headers: {
                Accept: "application/json",
            },
        });

        console.log("getAllUsers", res);
    };

    const addUser = async (username: string) => {
        // todo: Change the api endpoint
        const res = await fetch("add_user", {
            method: "POST",
            body: JSON.stringify({ username }),
        });

        console.log("addUser", res);
    };

    const updateSharedUsers = async (updatedUsers: any) => {
        // todo: Change the api endpoint
        const res = await fetch("update_users", {
            method: "POST",
            body: JSON.stringify({ users: updatedUsers }),
        });

        console.log("updateSharedUsers", res);
    };

    return {
        username,
        setUsername,
        sharedUsers,
        setSharedUsers,
        addUser,
        updateSharedUsers,
        getAllUsers,
    };
};

export default useProjectShare;
