import { describe, expect, test } from "vitest";
import { describeAnswer, diagnoseImportFailure, failedModuleUrl, probeUrl } from "@/utils/importFailure";

const URL_ = "https://ood.example/rnode/node/15476/flask/assets/BaseChart-CxpVgHMw.js";
const APACHE_502 =
    '<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML 2.0//EN">\n<html><head>\n<title>502 Proxy Error</title>\n</head><body>';

const answering = (body: string, init: ResponseInit) => async () => new Response(body, init);

describe("failedModuleUrl", () => {
    test.each([
        [`Failed to fetch dynamically imported module: ${URL_}`],
        [`error loading dynamically imported module: ${URL_}`],
        [`Unable to preload CSS for ${URL_}`],
    ])("finds the file in %s", (message) => {
        expect(failedModuleUrl(new TypeError(message))).toBe(URL_);
    });

    test("is undefined for errors that name no file", () => {
        expect(failedModuleUrl(new TypeError("Importing a module script failed."))).toBeUndefined();
        expect(failedModuleUrl(undefined)).toBeUndefined();
    });
});

describe("diagnoseImportFailure", () => {
    const error = new TypeError(`Failed to fetch dynamically imported module: ${URL_}`);

    test("reports a proxy's error page", async () => {
        const result = await diagnoseImportFailure(
            error,
            answering(APACHE_502, { status: 502, headers: { "content-type": "text/html; charset=iso-8859-1" } }),
        );
        expect(result?.probe).toMatchObject({
            url: URL_,
            status: 502,
            contentType: "text/html",
            title: "502 Proxy Error",
        });
        expect(result?.explanation).toContain("502 Proxy Error from a proxy");
    });

    test("says the file is there when fetching it again works", async () => {
        const result = await diagnoseImportFailure(
            error,
            answering("export {}", { status: 200, headers: { "content-type": "text/javascript" } }),
        );
        expect(result?.explanation).toMatch(/file is on the server/);
    });

    test("does not take a page served in place of the file for the file", async () => {
        const result = await diagnoseImportFailure(
            error,
            answering("<html><title>MDV</title></html>", { status: 200, headers: { "content-type": "text/html" } }),
        );
        expect(result?.explanation).toMatch(/got a web page \("MDV"\) instead of the file/);
    });

    test("points at a build mismatch for a missing file", async () => {
        const result = await diagnoseImportFailure(error, answering("nope", { status: 404 }));
        expect(result?.explanation).toMatch(/different build/);
    });

    test("reports a request that gets no answer", async () => {
        const failing = async (): Promise<Response> => {
            throw new TypeError("Failed to fetch");
        };
        const probe = await probeUrl(URL_, failing);
        expect(probe).toEqual({ url: URL_, fetchError: "Failed to fetch" });
    });

    test("is undefined when the error names no file", async () => {
        expect(await diagnoseImportFailure(new Error("something else"))).toBeUndefined();
    });
});

test("describeAnswer does not repeat a status the title already has", () => {
    expect(describeAnswer(502, "502 Proxy Error")).toBe("502 Proxy Error");
    expect(describeAnswer(500, "Internal Server Error")).toBe("500 (Internal Server Error)");
    expect(describeAnswer(404)).toBe("404");
});
