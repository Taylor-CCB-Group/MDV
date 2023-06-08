import DataStore from '../datastore/DataStore';
import { createEl } from '../utilities/Elements';

type ImageArrayConfig = {
    base_url: string,
    image_type: string,
    image_key: string,
    width: number,
    height: number,
    cancel?: boolean,
}

type ImageArrayEntry = {
    image: HTMLImageElement,
    zIndex: number,
}

/**
 * @class ImageArray 
 * @description
 * ImageArray is a class that holds a collection of images as a texture array, such that many images can be efficiently rendered
 * without requiring a lot of draw calls / GL state updates.
 * 
 */
export class ImageArray {
    textures: Map<string, ImageArrayEntry>;
    gl: WebGL2RenderingContext;
    logEl: HTMLDivElement;
    constructor(dataStore, canvas: HTMLCanvasElement, config: ImageArrayConfig) {
        this.textures = new Map();
        const gl = canvas.getContext("webgl2");
        this.gl = gl;
        this.logEl = createEl('div', {}, canvas.parentElement);
        this.logEl.style.color = "white";
        this.loadImageColumn(dataStore, gl, config);
    }

    drawProgress(n = 0) {
        this.logEl.textContent = `Loading images: ${Math.round(n * 100)}%`;
        this.gl.clearColor(n, n, n, 1);
        this.gl.clear(this.gl.COLOR_BUFFER_BIT);
    }
    loadImageColumn(ds: DataStore, gl: WebGL2RenderingContext, config: ImageArrayConfig) {
        // how about also passing a cancellation token or something?
        const { base_url, image_type, image_key, width, height } = config;
        const col = ds.columnIndex[image_key];
        const texture = gl.createTexture();
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D_ARRAY, texture);
        gl.texImage3D(gl.TEXTURE_2D_ARRAY, 0, gl.RGBA, col.data.length, width, height, 0, gl.RGBA, gl.UNSIGNED_BYTE, null);
        let nLoaded = 0;
        // this should be based on a selection filter (ie, we don't just call it on the whole column from constructor)
        col.data.map((d, i) => {
            //consider mapping from data index to texture index in case several rows have the same image?
            if (this.textures.has(d)) return;
            const url = `${base_url}/${d}.${image_type}`;
            const zIndex = i;
            const image = new Image(); //may want to use fetch instead for better control?
            this.textures.set(d, { image, zIndex });
            image.src = url;
            image.crossOrigin = "anonymous";
            image.onload = () => {
                if (config.cancel) return; // todo: cancel requests
                this.drawProgress(nLoaded++ / col.data.length);
                // bind texture and update it
                gl.activeTexture(gl.TEXTURE0);
                gl.bindTexture(gl.TEXTURE_2D_ARRAY, texture);
                //TODO: resize image to fit texture, release reference to image so it can be garbage collected
                gl.texSubImage3D(gl.TEXTURE_2D_ARRAY, 0, 0, 0, zIndex, image.width, image.height, 1, gl.RGBA, gl.UNSIGNED_BYTE, image);
                // revert to previous texture?
                // probably not going to hurt to leave it bound in this case
            }
        });
    }
}    
