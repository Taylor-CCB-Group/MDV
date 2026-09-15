import ProjectCard, { type ProjectCardProps } from "@/catalog/ProjectCard";
import { render, screen } from "@testing-library/react";
import { describe, expect, test, vi } from "vitest";

vi.mock("@/catalog/PermissionsContext", () => ({
    default: () => ({
        permissions: {
            createProject: true,
            importProject: true,
            deleteProject: true,
            renameProject: true,
            changeProjectAccess: true,
            exportProject: true,
            shareProject: true,
            editUserPermissions: true,
            removeUserFromProject: true,
        },
        isPublicPage: false,
    }),
}));

const commonProps: Omit<ProjectCardProps, "writable"> = {
    id: "10",
    name: "read-only-project",
    type: "Editable",
    lastModified: "2026-09-03 12:00:00",
    createdAt: "",
    owner: [],
    collaborators: [],
    numberOfStructures: "0",
    numberOfImages: "0",
    permissions: { read: true, edit: true, owner: true },
    onDelete: vi.fn(async () => {}),
    onRename: vi.fn(async () => {}),
    onChangeType: vi.fn(async () => {}),
    onAddCollaborator: vi.fn(),
    onExport: vi.fn(async () => {}),
};

describe("ProjectCard filesystem writability", () => {
    test("marks projects the server cannot write", () => {
        render(<ProjectCard {...commonProps} writable={false} />);

        expect(screen.getByText("Server read-only")).toBeDefined();
    });

    test("does not mark writable projects", () => {
        render(<ProjectCard {...commonProps} writable />);

        expect(screen.queryByText("Server read-only")).toBeNull();
    });
});
