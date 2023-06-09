import DataStore from '../datastore/DataStore';
import { DataModel } from '../table/DataModel';
import { createEl } from '../utilities/Elements';

export type ImageArrayConfig = {
    base_url: string,
    image_type: string,
    image_key: string,
    width: number,
    height: number,
    cancel?: boolean,
}

export type ImageArrayEntry = {
    // image: HTMLImageElement,
    zIndex: number,
    /** use in shader to scale uv coordinates */
    aspectRatio: number,
}

/**
 * @class ImageArray 
 * @description
 * ImageArray is a class that holds a collection of images as a texture array, such that many images can be efficiently rendered
 * without requiring a lot of draw calls / GL state updates.
 * 
 * Current design is brute-force, loading all images at once, little other logic.
 */
export class ImageArray {
    textures: Map<string, ImageArrayEntry>;
    texturesByIndex: Map<number, ImageArrayEntry>;
    texture: WebGLTexture;
    gl: WebGL2RenderingContext;
    logEl: HTMLDivElement;
    dataView: DataModel;
    onProgress: (n: number) => void;
    constructor(dataStore, canvas: HTMLCanvasElement, dataView: DataModel, config: ImageArrayConfig) {
        this.textures = new Map();
        this.texturesByIndex = new Map();
        
        this.dataView = dataView;

        const gl = canvas.getContext("webgl2");
        this.gl = gl;
        this.texture = gl.createTexture();
        this.logEl = createEl('div', {}, canvas.parentElement);
        this.logEl.style.color = "white";
        this.loadImageColumn(dataStore, gl, config);
    }
    getImageAspect(i: number) {
        return 0.5 + Math.random() // for testing
        return this.texturesByIndex.get(i).aspectRatio;
    }
    getImageIndex(i: number) {
        return this.texturesByIndex.get(i).zIndex;
    }
    drawProgress(n = 0) {
        this.logEl.textContent = `Loading images: ${Math.round(n * 100)}%`;
        this.onProgress && this.onProgress(n);
    }
    loadImageColumn(ds: DataStore, gl: WebGL2RenderingContext, config: ImageArrayConfig) {
        // how about also passing a cancellation token or something?
        const { base_url, image_type, image_key, width, height } = config;
        const col = ds.columnIndex[image_key];
        const texture = this.texture;
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D_ARRAY, texture);
        const mipLevels = Math.floor(Math.log2(width));

        gl.texParameteri(gl.TEXTURE_2D_ARRAY, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);
        gl.texParameteri(gl.TEXTURE_2D_ARRAY, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
        gl.texParameteri(gl.TEXTURE_2D_ARRAY, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
        gl.texParameteri(gl.TEXTURE_2D_ARRAY, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

        gl.texStorage3D(gl.TEXTURE_2D_ARRAY, mipLevels, gl.RGBA8, width, height, col.data.length);
        const memUsage = (width * height * 4 * col.data.length) / 1024 / 1024;
        //consider showing this in the UI ('i' for info?)
        console.log(`Allocated ${memUsage.toFixed(2)}MB for image array (not accounting for mipmaps)`);
        let nLoaded = 0;
        col.data.map((d, i: number) => {
            // XXX: still not working --- need to better grasp how DataModel works
            // this.dataView.updateModel();
            const imageName = d;//this.dataView.getItemField(d, image_key);
            if (this.textures.has(imageName)) {
                this.texturesByIndex.set(i, this.textures.get(imageName));
                nLoaded++;
                return;
            }
            const url = `${base_url}/${imageName}.${image_type}`;
            const zIndex = i;
            const image = new Image(); //may want to use fetch instead for better control?
            const entry = { zIndex, aspectRatio: 1 }; //not holding image ref, so we can garbage collect
            this.textures.set(imageName, entry);
            this.texturesByIndex.set(i, entry);
            image.src = url;
            image.crossOrigin = "anonymous";
            image.onload = () => {
                if (config.cancel) return; // todo: cancel requests
                // bind texture and update it
                gl.activeTexture(gl.TEXTURE0);
                gl.bindTexture(gl.TEXTURE_2D_ARRAY, texture);
                entry.aspectRatio = image.width / image.height; // mutating entry won't prompt re-computing vertex buffers...
                const resizedImage = resizeImage(image, width, height)//, this.logEl);
                gl.texSubImage3D(gl.TEXTURE_2D_ARRAY, 0, 0, 0, zIndex, width, height, 1, gl.RGBA, gl.UNSIGNED_BYTE, resizedImage);
                const error = gl.getError();
                if (error !== gl.NO_ERROR) {
                    console.error(`glError ${error} loading image #${zIndex} '${url}'`);
                }
                gl.generateMipmap(gl.TEXTURE_2D_ARRAY);
                // probably not going to hurt to leave it bound in this case
                this.drawProgress(nLoaded++ / col.data.length);
            }
            image.onerror = () => {
                console.error(`Error loading image #${zIndex} '${url}'`);
            }
        });
    }
}

function resizeImage(image: HTMLImageElement, width: number, height: number, logEl?: HTMLDivElement) {
    const canvas = document.createElement('canvas');
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