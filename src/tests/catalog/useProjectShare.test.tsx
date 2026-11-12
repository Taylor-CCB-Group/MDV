import { act, renderHook, waitFor } from "@testing-library/react";
import { beforeEach, describe, expect, test, vi } from "vitest";

const { apiFetchMock } = vi.hoisted(() => ({
    apiFetchMock: vi.fn(),
}));

vi.mock("@/utils/mdvRouting", () => ({
    apiFetch: apiFetchMock,
    buildApiUrl: (path: string) => path,
}));

import useProjectShare from "@/catalog/hooks/useProjectShare";

const ADMIN_ROLE_LOCKED =
    "An Admin is Owner of every project. Remove their Admin rights to change their project roles.";

function jsonResponse(body: unknown, status = 200) {
    return new Response(JSON.stringify(body), {
        status,
        headers: { "Content-Type": "application/json" },
    });
}

describe("useProjectShare refusals", () => {
    beforeEach(() => {
        apiFetchMock.mockReset();
        // The dialog loads its members when it opens.
        apiFetchMock.mockResolvedValueOnce(jsonResponse({ shared_users: [], all_users: [] }));
    });

    test("shows a refused removal next to the form and keeps the form", async () => {
        const { result } = renderHook(() => useProjectShare("10"));
        await waitFor(() => expect(apiFetchMock).toHaveBeenCalledTimes(1));
        apiFetchMock.mockResolvedValueOnce(jsonResponse({ error: ADMIN_ROLE_LOCKED }, 400));

        await act(async () => {
            await result.current.deleteSharedUser(1);
        });

        expect(result.current.errorMsg).toBe(ADMIN_ROLE_LOCKED);
        expect(result.current.error).toBe("");
    });

    test("still replaces the form when removing is forbidden", async () => {
        const { result } = renderHook(() => useProjectShare("10"));
        await waitFor(() => expect(apiFetchMock).toHaveBeenCalledTimes(1));
        apiFetchMock.mockResolvedValueOnce(jsonResponse({ error: "Only the project owner can remove users" }, 403));

        await act(async () => {
            await result.current.deleteSharedUser(1);
        });

        expect(result.current.error).toBe("Only the project owner can remove users");
    });
});
