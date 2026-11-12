import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";

const initializeFaro = vi.fn();

vi.mock("@grafana/faro-web-sdk", () => ({
    initializeFaro,
    getWebInstrumentations: () => [],
}));
vi.mock("@grafana/faro-web-tracing", () => ({
    TracingInstrumentation: class {},
}));

function serverSays(config: object) {
    const fetchMock = vi.fn().mockResolvedValue({ ok: true, json: async () => config });
    vi.stubGlobal("fetch", fetchMock);
    return fetchMock;
}

async function start() {
    // The module starts Faro at most once, so each test loads a fresh copy.
    const { startBrowserTelemetry } = await import("@/utils/browserTelemetry");
    await startBrowserTelemetry();
}

describe("startBrowserTelemetry", () => {
    beforeEach(() => {
        vi.resetModules();
        initializeFaro.mockClear();
        window.history.replaceState(null, "", "/project/1");
    });

    afterEach(() => {
        vi.unstubAllGlobals();
    });

    it("stays off when the server has no collector", async () => {
        const fetchMock = serverSays({ enabled: false });

        await start();

        expect(fetchMock).toHaveBeenCalledWith("/telemetry/config");
        expect(initializeFaro).not.toHaveBeenCalled();
    });

    it("posts to the deployment it was served from, under the server's name and sample rate", async () => {
        serverSays({ enabled: true, app_name: "mdv-test-browser", sample_rate: 0.5 });

        await start();

        expect(initializeFaro).toHaveBeenCalledTimes(1);
        const config = initializeFaro.mock.calls[0][0];
        expect(config.url).toBe(`${window.location.origin}/telemetry/faro`);
        expect(config.app.name).toBe("mdv-test-browser");
        expect(config.sessionTracking).toEqual({ samplingRate: 0.5 });
    });

    it("never starts when the API is on another origin", async () => {
        window.history.replaceState(null, "", "/?dir=https://elsewhere.example/project/1");
        const fetchMock = serverSays({ enabled: true, app_name: "x", sample_rate: 1 });

        await start();

        expect(fetchMock).not.toHaveBeenCalled();
        expect(initializeFaro).not.toHaveBeenCalled();
    });
});
