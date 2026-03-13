import { useEffect, useRef } from "react";
import { useChannelsStore, useViewerStore } from "./avivatorish/state";

type Range = [number, number];

type ThumbnailRenderer = {
    canvas: HTMLCanvasElement;
    gl: WebGL2RenderingContext;
    program: WebGLProgram;
    positionBuffer: WebGLBuffer;
    texture: WebGLTexture;
    positionLocation: number;
    textureLocation: WebGLUniformLocation;
    limitsLocation: WebGLUniformLocation;
    colorLocation: WebGLUniformLocation;
};

let sharedThumbnailRenderer: ThumbnailRenderer | null | undefined;

const clampRange = (range: Range, domain: Range): Range => [
    Math.max(domain[0], Math.min(range[0], domain[1])),
    Math.max(domain[0], Math.min(range[1], domain[1])),
];

const sortRange = ([start, end]: Range): Range => (start <= end ? [start, end] : [end, start]);

const createShader = (gl: WebGL2RenderingContext, type: number, source: string) => {
    const shader = gl.createShader(type);
    if (!shader) {
        throw new Error("Failed to create thumbnail shader");
    }
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        const info = gl.getShaderInfoLog(shader);
        gl.deleteShader(shader);
        throw new Error(info || "Failed to compile thumbnail shader");
    }
    return shader;
};

const createProgram = (gl: WebGL2RenderingContext, vertexSource: string, fragmentSource: string) => {
    const vertexShader = createShader(gl, gl.VERTEX_SHADER, vertexSource);
    const fragmentShader = createShader(gl, gl.FRAGMENT_SHADER, fragmentSource);
    const program = gl.createProgram();
    if (!program) {
        throw new Error("Failed to create thumbnail program");
    }
    gl.attachShader(program, vertexShader);
    gl.attachShader(program, fragmentShader);
    gl.linkProgram(program);
    gl.deleteShader(vertexShader);
    gl.deleteShader(fragmentShader);
    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
        const info = gl.getProgramInfoLog(program);
        gl.deleteProgram(program);
        throw new Error(info || "Failed to link thumbnail program");
    }
    return program;
};

const getSharedThumbnailRenderer = (): ThumbnailRenderer | null => {
    if (sharedThumbnailRenderer !== undefined) {
        return sharedThumbnailRenderer;
    }

    if (typeof document === "undefined") {
        sharedThumbnailRenderer = null;
        return sharedThumbnailRenderer;
    }

    try {
        const canvas = document.createElement("canvas");
        const gl = canvas.getContext("webgl2", {
            antialias: false,
            depth: false,
            premultipliedAlpha: false,
            preserveDrawingBuffer: true,
            stencil: false,
        });
        if (!gl) {
            sharedThumbnailRenderer = null;
            return sharedThumbnailRenderer;
        }

        const program = createProgram(
            gl,
            `#version 300 es
            in vec2 a_position;
            out vec2 v_uv;
            void main() {
                v_uv = a_position * 0.5 + 0.5;
                gl_Position = vec4(a_position, 0.0, 1.0);
            }`,
            `#version 300 es
            precision highp float;
            uniform sampler2D u_texture;
            uniform vec2 u_limits;
            uniform vec3 u_color;
            in vec2 v_uv;
            out vec4 outColor;
            void main() {
                float value = texture(u_texture, v_uv).r;
                float intensity = max(0.0, (value - u_limits.x) / max(0.0005, u_limits.y - u_limits.x));
                intensity = min(1.0, intensity);
                outColor = vec4(u_color * intensity, 1.0);
            }`,
        );

        const positionBuffer = gl.createBuffer();
        const texture = gl.createTexture();
        const textureLocation = gl.getUniformLocation(program, "u_texture");
        const limitsLocation = gl.getUniformLocation(program, "u_limits");
        const colorLocation = gl.getUniformLocation(program, "u_color");
        const positionLocation = gl.getAttribLocation(program, "a_position");

        if (
            !positionBuffer ||
            !texture ||
            !textureLocation ||
            !limitsLocation ||
            !colorLocation ||
            positionLocation < 0
        ) {
            sharedThumbnailRenderer = null;
            return sharedThumbnailRenderer;
        }

        gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);
        gl.bufferData(
            gl.ARRAY_BUFFER,
            new Float32Array([-1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, 1]),
            gl.STATIC_DRAW,
        );

        gl.bindTexture(gl.TEXTURE_2D, texture);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
        gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
        gl.pixelStorei(gl.PACK_ALIGNMENT, 1);

        sharedThumbnailRenderer = {
            canvas,
            gl,
            program,
            positionBuffer,
            texture,
            positionLocation,
            textureLocation,
            limitsLocation,
            colorLocation,
        };
        return sharedThumbnailRenderer;
    } catch (error) {
        console.error("Failed to initialize thumbnail renderer", error);
        sharedThumbnailRenderer = null;
        return sharedThumbnailRenderer;
    }
};

const normalizeSampleValue = (value: number, limits: Range) => {
    if (!Number.isFinite(value)) {
        return 0;
    }
    const [min, max] = limits;
    const scaled = Math.max(0, (value - min) / Math.max(0.0005, max - min));
    return Math.min(1, scaled);
};

const createFloatTextureData = (data: ArrayLike<number>, pixelCount: number) => {
    const values = new Float32Array(pixelCount);
    for (let i = 0; i < pixelCount; i += 1) {
        values[i] = Number(data[i] ?? 0);
    }
    return values;
};

const drawThumbnailWithCanvas2D = ({
    context,
    width,
    height,
    data,
    color,
    limits,
}: {
    context: CanvasRenderingContext2D;
    width: number;
    height: number;
    data: ArrayLike<number>;
    color: [number, number, number];
    limits: Range;
}) => {
    const pixelCount = Math.min(width * height, data.length);
    const imageData = context.createImageData(width, height);
    const pixels = imageData.data;

    for (let sourceIndex = 0; sourceIndex < pixelCount; sourceIndex += 1) {
        const pixelIndex = sourceIndex * 4;
        const intensity = normalizeSampleValue(Number(data[sourceIndex]), limits);

        pixels[pixelIndex] = Math.round(color[0] * intensity);
        pixels[pixelIndex + 1] = Math.round(color[1] * intensity);
        pixels[pixelIndex + 2] = Math.round(color[2] * intensity);
        pixels[pixelIndex + 3] = 255;
    }

    context.putImageData(imageData, 0, 0);
};

export default function ChannelThumbnail({ index }: { index: number }) {
    const raster = useChannelsStore((state) => state.raster[index]);
    const color = useChannelsStore((state) => state.colors[index] ?? [255, 255, 255]);
    const domain = useChannelsStore((state) => state.domains[index] ?? ([0, 1] as Range));
    const contrastLimits = useChannelsStore((state) => state.contrastLimits[index] ?? domain);
    const isChannelLoading = useViewerStore((state) => state.isChannelLoading[index]);
    const canvasRef = useRef<HTMLCanvasElement | null>(null);

    useEffect(() => {
        const canvas = canvasRef.current;
        if (!canvas || !raster?.width || !raster?.height || !raster.data) {
            return;
        }

        const { width, height, data } = raster;
        canvas.width = width;
        canvas.height = height;
        const clampedLimits = sortRange(clampRange(contrastLimits, domain));

        const displayContext = canvas.getContext("2d");
        if (!displayContext) {
            return;
        }

        displayContext.clearRect(0, 0, width, height);

        const renderer = getSharedThumbnailRenderer();
        if (!renderer) {
            drawThumbnailWithCanvas2D({
                context: displayContext,
                width,
                height,
                data,
                color,
                limits: clampedLimits,
            });
            return;
        }

        const pixelCount = Math.min(width * height, data.length);
        const textureData = createFloatTextureData(data, pixelCount);
        const { gl } = renderer;

        renderer.canvas.width = width;
        renderer.canvas.height = height;
        gl.viewport(0, 0, width, height);
        gl.useProgram(renderer.program);
        gl.bindBuffer(gl.ARRAY_BUFFER, renderer.positionBuffer);
        gl.enableVertexAttribArray(renderer.positionLocation);
        gl.vertexAttribPointer(renderer.positionLocation, 2, gl.FLOAT, false, 0, 0);

        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, renderer.texture);
        gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
        gl.pixelStorei(gl.PACK_ALIGNMENT, 1);
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.R32F, width, height, 0, gl.RED, gl.FLOAT, textureData);

        if (gl.getError() !== gl.NO_ERROR) {
            drawThumbnailWithCanvas2D({
                context: displayContext,
                width,
                height,
                data,
                color,
                limits: clampedLimits,
            });
            return;
        }

        gl.uniform1i(renderer.textureLocation, 0);
        gl.uniform2f(renderer.limitsLocation, clampedLimits[0], clampedLimits[1]);
        gl.uniform3f(renderer.colorLocation, color[0] / 255, color[1] / 255, color[2] / 255);
        gl.drawArrays(gl.TRIANGLES, 0, 6);

        displayContext.drawImage(renderer.canvas, 0, 0);
    }, [color, contrastLimits, domain, raster]);

    if (isChannelLoading || !raster?.width || !raster?.height || !raster.data?.length) {
        return (
            <div className="flex h-[68px] w-full items-center justify-center rounded-md border border-dashed border-[hsl(var(--border))] bg-[hsl(var(--muted))] px-2 text-[10px] font-medium uppercase tracking-[0.14em] text-[hsl(var(--muted-foreground))]">
                Sample unavailable
            </div>
        );
    }

    const aspectRatio = raster.width / raster.height;

    return (
        <div className="flex w-full items-center justify-center overflow-hidden rounded-md border border-[hsl(var(--border))] bg-[radial-gradient(circle_at_top,hsl(var(--muted)/0.35),hsl(var(--background)))] p-1">
            <canvas
                ref={canvasRef}
                className="block h-auto w-full max-w-full rounded-sm"
                style={{
                    aspectRatio: `${raster.width} / ${raster.height}`,
                    imageRendering: "pixelated",
                    maxHeight: aspectRatio > 3 ? "88px" : "120px",
                    objectFit: "contain",
                }}
            />
        </div>
    );
}
