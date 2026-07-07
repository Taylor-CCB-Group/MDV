import type { Dispatch, SetStateAction } from "react";
import { FolderKanban, UserPlus, UserRound, type LucideIcon } from "lucide-react";
import type { AdminPermission, AdminUser } from "./api";

export type InitialProjectAccess = {
    projectId: number;
    permission: AdminPermission;
};

export type AdminAccessState = "login_required" | "admin_required" | null;
export type AdminRoute = "users" | "projects" | "profile";

export type PermissionSetter = Dispatch<SetStateAction<AdminPermission>>;

export const adminRoutes: Array<{ id: AdminRoute; label: string; description: string; icon: LucideIcon }> = [
    {
        id: "users",
        label: "Manage Users",
        description: "Create deployment users and review account status",
        icon: UserPlus,
    },
    {
        id: "projects",
        label: "Manage Projects",
        description: "Manage project membership and permission levels",
        icon: FolderKanban,
    },
    {
        id: "profile",
        label: "Profile",
        description: "Review the current admin session",
        icon: UserRound,
    },
];

export function formatUserName(user: AdminUser) {
    const name = [user.firstName, user.lastName].filter(Boolean).join(" ");
    return name || user.email;
}

export function getRouteDescription(route: AdminRoute) {
    return adminRoutes.find((item) => item.id === route)?.description ?? adminRoutes[0].description;
}

export function permissionFromValue(value: string): AdminPermission | null {
    if (value === "view" || value === "edit" || value === "owner") {
        return value;
    }
    return null;
}
