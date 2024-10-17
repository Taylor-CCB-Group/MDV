import type { Device, Texture } from "@luma.gl/core";
import { getProjectURL } from "../dataloaders/DataLoaderUtil";
import type DataStore from "../datastore/DataStore";
import type { DataModel } from "../table/DataModel";
import { createEl } from "../utilities/Elements";

export type ImageArrayConfig = {
    base_url: string;
    image_type: string;
    image_key: string;
    width: number;
    height: number;
    cancel?: boolean;
};

export type ImageArrayEntry = {
    // image: HTMLImageElement,
    zIndex: number;
    /** use in shader to scale uv coordinates */
    aspectRatio: number;
    url: string;
};

/**
 * @class ImageArray
 * @description
 * ImageArray is a class that holds a collection of images as a texture array, such that many images can be efficiently rendered
 * without requiring a lot of draw calls / GL state updates.
 *
 * Current design is brute-force, loading all images at once, little other logic.
 * Will fall over and die if you try to use more than "Max Array Texture Layers" (usually 2048) images.
 */
export class ImageArray {
    textures: Map<string, ImageArrayEntry>;
    texturesByIndex: Map<number, ImageArrayEntry>;
    texture: WebGLTexture;
    gl: WebGL2RenderingContext;
    logEl: HTMLElement;
    dataView: DataModel;
    lumaTexture?: Texture;
    onProgress: (n: number) => void;
    constructor(
        dataStore,
        canvas: HTMLCanvasElement,
        dataView: DataModel,
        config: ImageArrayConfig,
    ) {
        this.textures = new Map();
        this.texturesByIndex = new Map();

        this.dataView = dataView;

        const gl = canvas.getContext("webgl2");
        this.gl = gl;
        this.texture = gl.createTexture();
        this.logEl = createEl("div", {}, canvas.parentElement);
        this.logEl.style.color = "white";
        this.loadImageColumn(dataStore, gl, config);
    }
    wrapLumaTexture(device: Device) {
        if (this.lumaTexture) return;
        console.warn("wrapping ImageArray texture as luma texture - todo: refactor in future");
        this.lumaTexture = device.createTexture({
            dimension: '2d-array',
        });
        (this.lumaTexture as any).handle = this.texture;
    }
    getImageAspect(i: number) {
        // return 0.5 + Math.random() // for testing
        return this.texturesByIndex.get(i).aspectRatio;
    }
    getImageIndex(i: number) {
        return this.texturesByIndex.get(i).zIndex;
    }
    drawProgress(n = 0) {
        // this.logEl.textContent = `Loading images: ${Math.round(n * 100)}%`;
        this.onProgress?.(n);
    }
    loadImageColumn(
        ds: DataStore,
        gl: WebGL2RenderingContext,
        config: ImageArrayConfig,
    ) {
        // how about also passing a cancellation token or something?
        const { image_type, image_key, width, height } = config;
        const base_url = getProjectURL(config.base_url);
        const col = ds.columnIndex[image_key];
        const texture = this.texture;
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D_ARRAY, texture);
        const mipLevels = Math.floor(Math.log2(width));

        gl.texParameteri(
            gl.TEXTURE_2D_ARRAY,
            gl.TEXTURE_MIN_FILTER,
            gl.LINEAR_MIPMAP_LINEAR,
        );
        gl.texParameteri(gl.TEXTURE_2D_ARRAY, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
        gl.texParameteri(
            gl.TEXTURE_2D_ARRAY,
            gl.TEXTURE_WRAP_S,
            gl.CLAMP_TO_EDGE,
        );
        gl.texParameteri(
            gl.TEXTURE_2D_ARRAY,
            gl.TEXTURE_WRAP_T,
            gl.CLAMP_TO_EDGE,
        );

        const isUnique = col.datatype === "unique";
        const numImages: number = isUnique
            ? col.data.length / col.stringLength
            : col.values?.length || col.data.length;

        gl.texStorage3D(
            gl.TEXTURE_2D_ARRAY,
            mipLevels,
            gl.RGBA8,
            width,
            height,
            numImages,
        );
        const memUsage = (width * height * 4 * numImages) / 1024 / 1024;
        //consider showing this in the UI ('i' for info?)
        console.log(
            `Allocated ${memUsage.toFixed(2)}MB for image array (not accounting for mipmaps)`,
        );
        let nLoaded = 0;
        function getTextArray(): string[] {
            if (!isUnique) {
                if (col.values === undefined) return col.data; //not really a string array, but whatever
                return [...col.data].map((d) => col.values[d]); //slightly garbagey, never mind
            }
            const tc = new TextDecoder();
            const textArray = new Array(numImages);
            for (let i = 0; i < numImages; i++) {
                const start = i * col.stringLength;
                const end = start + col.stringLength;
                textArray[i] = tc.decode(col.data.slice(start, end));
            }
            return textArray;
        }
        const data = getTextArray();
        data.forEach((d, i: number) => {
            // XXX: still not working --- need to better grasp how DataModel works
            // this.dataView.updateModel();
            // -> would be better to map over col.values in the first place?
            const imageName = d; //col.values[d];//this.dataView.getItemField(d, image_key);
            if (this.textures.has(imageName)) {
                this.texturesByIndex.set(i, this.textures.get(imageName));
                nLoaded++;
                return;
            }
            const url = `${base_url}/${imageName}.${image_type}`;
            const zIndex = i;
            const image = new Image(); //may want to use fetch instead for better control?
            const entry = { zIndex, aspectRatio: 1, url }; //not holding image ref, so we can garbage collect
            this.textures.set(imageName, entry);
            this.texturesByIndex.set(i, entry);
            image.src = url;
            image.crossOrigin = "anonymous";
            image.onload = () => {
                if (config.cancel) return; // todo: cancel requests
                // bind texture and update it
                gl.activeTexture(gl.TEXTURE0);
                gl.bindTexture(gl.TEXTURE_2D_ARRAY, texture);
                // TODO consider changing how aspect ratio works to allow for non-square texture
                entry.aspectRatio = image.width / image.height; // mutating entry won't prompt re-computing vertex buffers...
                const resizedImage = resizeImage(image, width, height); //, this.logEl);
                gl.texSubImage3D(
                    gl.TEXTURE_2D_ARRAY,
                    0,
                    0,
                    0,
                    zIndex,
                    width,
                    height,
                    1,
                    gl.RGBA,
                    gl.UNSIGNED_BYTE,
                    resizedImage,
                );
                const error = gl.getError();
                if (error !== gl.NO_ERROR) {
                    console.error(
                        `glError ${error} loading image #${zIndex} '${url}'`,
                    );
                }
                gl.generateMipmap(gl.TEXTURE_2D_ARRAY);
                this.drawProgress(nLoaded++ / numImages);
            };
            image.onerror = () => {
                console.error(`Error loading image #${zIndex} '${url}'`);
            };
        });
    }
}

function resizeImage(
    image: HTMLImageElement,
    width: number,
    height: number,
    logEl?: HTMLDivElement,
) {
    const canvas = document.createElement("canvas");
    canvas.width = width;
    canvas.height = height;
    const ctx = canvas.getContext("2d");
    ctx.drawImage(image, 0, 0, width, height);
    if (logEl) {
        canvas.style.opacity = "0.2";
        logEl.parentElement.appendChild(canvas);
        setTimeout(() => {
            logEl.parentElement.removeChild(canvas);
        }, 1000);
    }
    return canvas;
}
