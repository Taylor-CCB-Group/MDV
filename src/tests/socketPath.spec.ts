import { beforeEach, describe, expect, test } from "vitest";
import { getSocketPath } from "@/utils/file_upload/fileUploadHelpers";

function setUrl(path: string) {
    window.history.replaceState({}, "", path);
}

describe("getSocketPath", () => {
    beforeEach(() => {
        setUrl("/");
    });

    // state.json sends mdv_api_root "/", meaning "wherever the app is mounted". A proxy
    // that strips a path prefix, such as Verily Workbench's, is invisible to the server,
    // so only the page location knows the mount path.
    test("keeps the mount path of a project page behind a prefix-stripping proxy", () => {
        setUrl("/app/abc/proxy/5055/project/1");

        expect(getSocketPath("/")).toBe("/app/abc/proxy/5055/socket.io");
    });

    test("keeps the mount path when the project page has a trailing slash", () => {
        setUrl("/app/abc/proxy/5055/project/1/");

        expect(getSocketPath("/")).toBe("/app/abc/proxy/5055/socket.io");
    });

    test("uses the root for an app served at the root", () => {
        setUrl("/project/1");

        expect(getSocketPath("/")).toBe("/socket.io");
    });

    test("uses an explicit api root", () => {
        setUrl("/project/1");

        expect(getSocketPath("/mdv/")).toBe("/mdv/socket.io");
    });

    test("takes only the path from an absolute api root", () => {
        setUrl("/project/1");

        expect(getSocketPath("http://localhost:5055/mdv/")).toBe("/mdv/socket.io");
    });
});
