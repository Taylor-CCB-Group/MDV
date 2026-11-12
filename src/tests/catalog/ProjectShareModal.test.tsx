import { fireEvent, render, screen, within } from "@testing-library/react";
import { beforeEach, describe, expect, test, vi } from "vitest";
import type { RegisteredUser, SharedUser } from "@/catalog/hooks/useProjectShare";
import ProjectShareModal from "@/catalog/ProjectShareModal";

// An Admin is Owner of every project, so the server refuses to change or
// remove their role, and adds them only as Owner (ADR 0019).

vi.mock("@/catalog/PermissionsContext", () => ({
    default: () => ({
        permissions: { editUserPermissions: true, removeUserFromProject: true },
        isPublicPage: false,
    }),
}));

const addUser = vi.fn(async () => {});
let sharedUsers: SharedUser[] = [];
let newUser: RegisteredUser | null = null;

vi.mock("@/catalog/hooks/useProjectShare", () => ({
    default: () => ({
        email: "",
        setEmail: vi.fn(),
        sharedUsers,
        addUser,
        isLoading: false,
        error: "",
        errorMsg: "",
        userList: [],
        newUser,
        setNewUser: vi.fn(),
        deleteSharedUser: vi.fn(async () => {}),
        changeUserPermission: vi.fn(async () => {}),
    }),
}));

function rowFor(email: string) {
    const row = screen.getByText(email).closest("tr");
    if (!row) throw new Error(`no row for ${email}`);
    return within(row);
}

describe("ProjectShareModal and Admins", () => {
    beforeEach(() => {
        addUser.mockClear();
        sharedUsers = [
            { id: 1, email: "admin@example.org", permission: "Owner", is_admin: true },
            { id: 2, email: "ann@example.org", permission: "View", is_admin: false },
        ];
        newUser = null;
    });

    test("locks an Admin's role and leaves other members' roles editable", () => {
        render(<ProjectShareModal open onClose={vi.fn()} projectId="10" />);

        const admin = rowFor("admin@example.org");
        expect(admin.getByText("Admin")).toBeDefined();
        expect(admin.getByRole("combobox").getAttribute("aria-disabled")).toBe("true");
        expect(admin.getByRole("button", { name: "Remove" }).hasAttribute("disabled")).toBe(true);

        const member = rowFor("ann@example.org");
        expect(member.getByRole("combobox").getAttribute("aria-disabled")).toBeNull();
        expect(member.getByRole("button", { name: "Remove" }).hasAttribute("disabled")).toBe(false);
    });

    test("adds an Admin as Owner and anyone else as View", () => {
        newUser = { id: 3, email: "other-admin@example.org", is_admin: true };
        const { rerender } = render(<ProjectShareModal open onClose={vi.fn()} projectId="10" />);
        fireEvent.click(screen.getByRole("button", { name: "Add" }));
        expect(addUser).toHaveBeenLastCalledWith(3, "Owner");

        newUser = { id: 4, email: "bo@example.org", is_admin: false };
        rerender(<ProjectShareModal open onClose={vi.fn()} projectId="10" />);
        fireEvent.click(screen.getByRole("button", { name: "Add" }));
        expect(addUser).toHaveBeenLastCalledWith(4, "View");
    });
});
