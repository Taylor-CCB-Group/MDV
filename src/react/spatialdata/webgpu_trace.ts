import { webgpuTraceSetting } from "@/lib/browserSettings";

/**
 * Opt-in tracing of WebGPU command encoder use, for errors whose cause the error itself
 * cannot show: validation errors arrive asynchronously, so their stack is the reporting
 * path, not the caller. ("GPUCommandEncoder.finish: encoder state is 'Locked'" on Safari
 * was the first.) Turned on by "WebGPU call tracing" in the Debug dialog's browser
 * settings, from the next page load: it has to wrap the API before any device exists.
 *
 * It wraps the WebGPU prototypes themselves, so it sees every caller (luma, intraspatial's
 * raw encoders), and records:
 * - a pass still open on an encoder when it is finished or another pass begins on it,
 *   with both stacks;
 * - for each uncaptured device error, the encoder calls just before it. If none of them is
 *   at fault, the browser's own encoding is the likely suspect.
 *
 * Findings go to the console and to a chart's "Copy details" (`webgpuTraceReport`).
 */

type TracedCall = { call: string; label: string; stack: string };

const RECENT_CALLS = 40;
const MAX_FINDINGS = 20;
const recent: TracedCall[] = [];
const findings: string[] = [];

/**
 * The caller's frames, without the "Error" header Chrome adds (Safari has none) and
 * without this module's own: `depth` is how many of them sit above the caller, counting
 * this function.
 */
function callerStack(depth: number): string {
    return (new Error().stack ?? "")
        .split("\n")
        .filter((line) => line.trim() && !/^Error\b/.test(line))
        .slice(depth, depth + 5)
        .map((line) => line.trim())
        .join("\n");
}

function labelOf(object: { label: string }) {
    return object.label ? `"${object.label}"` : "(unlabeled)";
}

function record(call: string, label: string, stack: string) {
    recent.push({ call, label, stack });
    if (recent.length > RECENT_CALLS) recent.shift();
}

function report(finding: string) {
    console.error(`[MDV WebGPU trace] ${finding}`);
    if (findings.length < MAX_FINDINGS) findings.push(finding);
}

function formatRecent(count = RECENT_CALLS) {
    return recent
        .slice(-count)
        .map(({ call, label, stack }) => `${call} ${label}\n${stack.replace(/^/gm, "    ")}`)
        .join("\n");
}

function install() {
    const open = new WeakMap<GPUCommandEncoder, string>();

    const beginning = (encoder: GPUCommandEncoder, call: string) => {
        const stack = callerStack(3);
        const begunAt = open.get(encoder);
        if (begunAt) {
            report(`${call} on encoder ${labelOf(encoder)} while a pass on it is open.\nPass begun at:\n${begunAt}\n${call} at:\n${stack}`);
        }
        record(call, labelOf(encoder), stack);
        open.set(encoder, stack);
    };
    const trackEnd = (encoder: GPUCommandEncoder, pass: GPURenderPassEncoder | GPUComputePassEncoder) => {
        const end = pass.end;
        pass.end = () => {
            open.delete(encoder);
            record("pass.end", labelOf(encoder), callerStack(2));
            end.call(pass);
        };
    };

    const encoderProto = GPUCommandEncoder.prototype;
    const beginRenderPass = encoderProto.beginRenderPass;
    encoderProto.beginRenderPass = function (this: GPUCommandEncoder, descriptor: GPURenderPassDescriptor) {
        beginning(this, "beginRenderPass");
        const pass = beginRenderPass.call(this, descriptor);
        trackEnd(this, pass);
        return pass;
    };
    const beginComputePass = encoderProto.beginComputePass;
    encoderProto.beginComputePass = function (this: GPUCommandEncoder, descriptor?: GPUComputePassDescriptor) {
        beginning(this, "beginComputePass");
        const pass = beginComputePass.call(this, descriptor);
        trackEnd(this, pass);
        return pass;
    };
    const finish = encoderProto.finish;
    encoderProto.finish = function (this: GPUCommandEncoder, descriptor?: GPUCommandBufferDescriptor) {
        const stack = callerStack(2);
        const begunAt = open.get(this);
        if (begunAt) {
            report(`finish() on encoder ${labelOf(this)} while a pass on it is open.\nPass begun at:\n${begunAt}\nfinish() at:\n${stack}`);
        }
        open.delete(this);
        record("finish", labelOf(this), stack);
        return finish.call(this, descriptor);
    };

    const submit = GPUQueue.prototype.submit;
    GPUQueue.prototype.submit = function (this: GPUQueue, commandBuffers: Iterable<GPUCommandBuffer>) {
        record("queue.submit", labelOf(this), callerStack(2));
        return submit.call(this, commandBuffers);
    };

    const requestDevice = GPUAdapter.prototype.requestDevice;
    GPUAdapter.prototype.requestDevice = async function (this: GPUAdapter, descriptor?: GPUDeviceDescriptor) {
        const device = await requestDevice.call(this, descriptor);
        device.addEventListener("uncapturederror", (event) => {
            report(`Uncaptured error on device ${labelOf(device)}: ${event.error.message}\nRecent encoder calls, oldest first:\n${formatRecent(8)}`);
        });
        return device;
    };

    console.warn("[MDV] WebGPU call tracing is on (Debug dialog > Browser settings). It slows rendering a little.");
}

export const webgpuTraceEnabled = webgpuTraceSetting.get() === "on" && typeof GPUCommandEncoder !== "undefined";
if (webgpuTraceEnabled) install();

/** What tracing found, for a chart's "Copy details"; null while tracing is off. */
export function webgpuTraceReport(): string | null {
    if (!webgpuTraceEnabled) return null;
    if (findings.length) return findings.join("\n\n");
    return `Nothing found. Last encoder calls, oldest first:\n${formatRecent(8)}`;
}
