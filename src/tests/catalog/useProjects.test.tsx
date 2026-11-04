import { act, renderHook } from "@testing-library/react";
import { beforeEach, describe, expect, test, vi } from "vitest";

const { apiFetchMock } = vi.hoisted(() => ({
    apiFetchMock: vi.fn(),
}));

vi.mock("@/utils/mdvRouting", () => ({
    apiFetch: apiFetchMock,
    buildApiUrl: (path: string) => path,
}));

import useProjects from "@/catalog/hooks/useProjects";

function jsonResponse(body: unknown) {
    return new Response(JSON.stringify(body), {
        status: 200,
        headers: { "Content-Type": "application/json" },
    });
}

describe("useProjects writability", () => {
    beforeEach(() => {
        apiFetchMock.mockReset();
    });

    test("maps filesystem writability from the projects response", async () => {
        apiFetchMock.mockResolvedValueOnce(
            jsonResponse([
                {
                    id: 10,
                    name: "read-only-project",
                    writable: false,
                },
            ]),
        );
        const { result } = renderHook(() => useProjects());

        await act(async () => {
            await result.current.fetchProjects();
        });

        expect(result.current.projects[0].writable).toBe(false);
    });

    test("reports projects registered without write access", async () => {
        apiFetchMock
            .mockResolvedValueOnce(
                jsonResponse({
                    created_project_ids: [10],
                    unwritable_projects: [{ id: 10, name: "read-only-project" }],
                }),
            )
            .mockResolvedValueOnce(jsonResponse([]));
        const { result } = renderHook(() => useProjects());

        await act(async () => {
            await result.current.rescanProjects();
        });

        expect(result.current.rescanWarning).toContain("read-only-project");
        expect(result.current.rescanWarning).toContain("will open read-only");
    });
});

describe("useProjects restore", () => {
    beforeEach(() => {
        apiFetchMock.mockReset();
    });

    test("fetches only the restored project, not the whole list", async () => {
        apiFetchMock
            .mockResolvedValueOnce(jsonResponse([{ id: 1, name: "kept" }]))
            .mockResolvedValueOnce(jsonResponse({ restoredProjectId: 5 }))
            .mockResolvedValueOnce(jsonResponse({ id: 5, name: "restored", writable: true }));
        const { result } = renderHook(() => useProjects());

        await act(async () => {
            await result.current.fetchProjects();
        });
        await act(async () => {
            await result.current.restoreRecycledProject("5");
        });

        expect(apiFetchMock.mock.calls.map(([path]) => path)).toEqual([
            "projects",
            "projects/recycle-bin/5/restore",
            "projects/5",
        ]);
        expect(result.current.projects.map((p) => p.name).sort()).toEqual(["kept", "restored"]);
    });

    test("falls back to the whole list if the restored project can't be fetched", async () => {
        apiFetchMock
            .mockResolvedValueOnce(jsonResponse({ restoredProjectId: 5 }))
            .mockResolvedValueOnce(new Response("{}", { status: 404 }))
            .mockResolvedValueOnce(jsonResponse([{ id: 5, name: "restored" }]));
        const { result } = renderHook(() => useProjects());

        await act(async () => {
            await result.current.restoreRecycledProject("5");
        });

        expect(apiFetchMock.mock.calls.map(([path]) => path)).toEqual([
            "projects/recycle-bin/5/restore",
            "projects/5",
            "projects",
        ]);
        expect(result.current.projects.map((p) => p.name)).toEqual(["restored"]);
    });

    test("falls back to the whole list if fetching the restored project throws", async () => {
        const consoleError = vi.spyOn(console, "error").mockImplementation(() => {});
        apiFetchMock
            .mockResolvedValueOnce(jsonResponse({ restoredProjectId: 5 }))
            .mockRejectedValueOnce(new TypeError("Failed to fetch"))
            .mockResolvedValueOnce(jsonResponse([{ id: 5, name: "restored" }]));
        const { result } = renderHook(() => useProjects());

        try {
            await act(async () => {
                await result.current.restoreRecycledProject("5");
            });
        } finally {
            consoleError.mockRestore();
        }

        expect(result.current.isErrorModalOpen).toBe(false);
        expect(result.current.projects.map((p) => p.name)).toEqual(["restored"]);
    });
});
