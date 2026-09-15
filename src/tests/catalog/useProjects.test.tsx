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
