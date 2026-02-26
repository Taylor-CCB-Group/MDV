import regl from "regl";
import * as glMatrix from "gl-matrix";
import { createEl, makeDraggable } from "../utilities/Elements.js";
import { Camera } from "./Camera.js";

class WGL2DI {
    /**
     * Creates a new wg2di instance
     * @param {string|object} div - The id of the div to house the instance or jquery element
     * @param {object} config (optional)
     * - in_view_only - if true then only those objects in view will be drawn after panning/zooming.This will
     * speed things up if maniluplating individual objects but slow down panning and zooming
     * - allow_object_drag - if true individual objects can be dragged 
     * - default_brush - if true then brushing will be the default drag behavior - uses shift for panning
     *  Otherwise panning is default and shift used for brushing
     * - noCameraControl - if true then it won't use pan/zoom controls
     */
    constructor(div, config) {
        if (!config) {
            config = {};
        }
        this.config = config;
        if (!config.offset) {
            config.offset = [0, 0, 0, 0];
        }
        this.mode = config.mode || "2d";
        this.movingImage = 0;
        this.__doc__ = document;

        this.config.brush = this.config.brush || "default";

        /*	this.draw_options=config.draw_options?config.draw_options:{depth:{enable:true},
		blend:{}
	}
	*/

        this.draw_options = config.draw_options
            ? config.draw_options
            : {
                  depth: { enable: false },
                  blend: {
                      enable: true,
                      func: {
                          srcRGB: "src alpha",
                          srcAlpha: "src alpha",
                          dstRGB: "one minus src alpha",
                          dstAlpha: "one minus src alpha",
                      },
                  },
              };

        this.regl = null;
        this.pickbuffer = null;

        this.hideOnFilter = true;

        this.pointOpacity = 0.8;
        this.pointRadius = 10;
        this.isFiltered = false;

        this.x_scale = 1.0;
        this.y_scale = 1.0;
        this.offset = [0, 0];

        //html elements
        if (typeof div === "string") {
            this.div_container = document.getElementById(div);
        } else {
            this.div_container = div;
        }

        this._setUpDocument(config.width, config.height);

        //handlers
        this.handlers = {
            object_clicked: {},
            object_over: {},
            object_out: {},
            brush_stopped: {},
            zoom_stopped: {},
            panning_stopped: {},
            pan_or_zoom: {},
        };

        //switches
        this.draw_labels = false;

        this.object_types = [];
        this.pointScale = 0;

        // circle shapes
        this.circle_properties = {
            x_pos: 1,
            y_pos: 1,
            color: 1,
            pick_color: 1,
            localFilter: 1,
            globalFilter: 1,
            colorFilter: 1,
        };

        this.circles = {};
        for (const prop in this.circle_properties) {
            this.circles[prop] = [];
        }
        this.circles.count = 0;
        this.object_types.push({
            data: this.circles,
            properties: this.circle_properties,
            vertices: 1,
            primitive: "points",
        });

        // line shapes
        this.line_properties = { position: 2, color: 2, opacity: 2 };

        this.lines = {};
        for (const prop in this.line_properties) {
            this.lines[prop] = [];
        }
        this.lines.count = 0;
        this.object_types.push({
            data: this.lines,
            properties: this.line_properties,
            vertices: 2,
            primitive: "lines",
        });

        //images
        this.image_properties = { position: 2, color: 2, opacity: 2 };
        this.images = {};
        for (const prop in this.image_properties) {
            this.images[prop] = [];
        }
        this.images.count = 0;

        //special propertis for images
        this.images.props = { x_y: [], w_h: [], text: [] };

        this.object_types.push({
            data: this.images,
            properties: this.image_properties,
            vertices: 6,
            primitive: "triangles",
        });

        if (this.mode === "3d") {
            this.circle_properties.z_pos = 1;
            this.camera = new Camera({
                center: config.cameraCenter,
                distance: config.cameraDistance || 500,
                theta: 0.75,
                phi: 0.5,
            });
            this.axisScales = [1, 1, 1];
        }

        //squares
        this.square_properties = {
            x_pos: 1,
            y_pos: 1,
            color: 1,
            pick_color: 1,
            localFilter: 1,
            globalFilter: 1,
        };

        this.squares = {};
        for (const prop in this.square_properties) {
            this.squares[prop] = [];
        }
        this.squares.count = 0;
        this.object_types.push({
            data: this.squares,
            properties: this.square_properties,
            vertices: 1,
            primitive: "points",
        });

        this.rectangle_properties = { position: 6, color: 6, opacity: 6 };

        this.rects = {};
        for (const prop in this.rectangle_properties) {
            this.rects[prop] = [];
        }
        this.rects.count = 0;
        this.object_types.push({
            data: this.rects,
            properties: this.rectangle_properties,
            vertices: 6,
            primitive: "triangles",
        });

        //The last mouse position recorded
        this.mouse_position = null;
        //Was an object clicked
        this.object_clicked = null;
        //an object was clicked
        this.dragging = false;
        //object which mouse is over
        this.object_mouse_over = null;

        this.zoom_amount = 0;
        this.draw_order =
            this.mode === "3d" ? [2, 0, 1, 3, 4] : [2, 0, 1, 3, 4];

        regl({
            onDone: (err, regl) => {
                this.regl = regl;
                regl._refresh();
                this.pickbuffer = regl.framebuffer({
                    colorFormat: "rgba",
                    height: this.height,
                    width: this.width,
                });
                this._initDrawMethods();
                if (this.mode === "3d") {
                    this.object_types[0]["method"] = this.__draw3DCircles;
                    this.object_types[1]["method"] = this.__draw3DLines;
                }

                this._addHandlers();
            },
            canvas: this.canvas,
            attributes: {
                antialias: false,
                preserveDrawingBuffer: true, //allows us to screenshot entire page without app.refresh() logic.
                willReadFrequently: true,
            },
        });
    }

    setBackGroundColor(color) {
        color = color || "";
        this.div_container.style.background = color;
    }

    setCameraProperty(property, value) {
        this.camera.cameraState[property] = value;
        this.camera.updateCamera();
    }

    _setUpDocument(width, height) {
        if (this.config.brush) {
            this.div_container.style.cursor = "crosshair";
        }
        if (!height) {
            const box = this.div_container.getBoundingClientRect();
            this.height = height = Math.round(box.height);
            this.width = width = Math.round(box.width);
        } else {
            this.height = height;
            this.width = width;
            this.div_container.style.height = `${height}px`;
            this.div_container.style.width = `${width}px`;
        }

        const attr = {
            height: height,
            width: width,
            styles: {
                position: "absolute",
                left: "0px",
                top: "0px",
            },
        };
        this.canvas = createEl("canvas", attr);
        this.label_canvas = createEl("canvas", attr);
        this.label_context = this.label_canvas.getContext("2d");
        this.div_container.append(this.canvas);
        this.div_container.append(this.label_canvas);
    }

    colorPoints(colorFunc) {
        const len = this.circles.count;
        const colors = this.circles.color;
        const cf = this.circles.colorFilter;
        for (let n = 0; n < len; n++) {
            const col = colorFunc(n);
            const cp = n * 3;
            if (col) {
                colors[cp] = col[0];
                colors[cp + 1] = col[1];
                colors[cp + 2] = col[2];
                cf[n] = 1;
            } else {
                cf[n] = 0;
            }
        }
        //color highlights
        if (this.highlightPoints) {
            const ps = this.highlightPoints;
            for (let i = 0; i < ps.length; i++) {
                const ci = ps.indexes[i] * 3;
                const hi = i * 3;
                ps.color[hi] = this.circles.color[ci];
                ps.color[hi + 1] = this.circles.color[ci + 1];
                ps.color[hi + 2] = this.circles.color[ci + 2];
            }
        }
    }

    remove() {
        this.pickbuffer.destroy();
        this.regl.destroy();
        //this.canvas.attr({height:1,width:1});
        this.canvas.remove();
        this.label_canvas.remove();
    }

    setSize(width, height, rescale = null) {
        if (this.height === height && this.width === width) {
            return;
        }
        width = Math.round(width);
        height = Math.round(height);
        if (rescale) {
            const scale = rescale(width, height);
            this.x_scale = scale[0] || this.x_scale;
            this.y_scale = scale[1] || thisy_scale;
        } else {
            const x_ratio = width / this.width;
            const y_ratio = height / this.height;

            this.x_scale = this.x_scale * x_ratio;
            this.y_scale = this.y_scale * y_ratio;
        }

        this.height = height;
        this.width = width;
        this.div_container.style.height = `${height}px`;
        this.div_container.style.width = `${width}px`;
        this.canvas.setAttribute("height", height);
        this.canvas.setAttribute("width", width);
        this.label_canvas.setAttribute("height", height);
        this.label_canvas.setAttribute("width", width);

        try {
            this.pickbuffer.destroy();
        } catch (e) {
            // we mostly avoid situations like this, which was arising with -ve width/height leading to inconsistent state.
            console.error(
                `caught exception '${e}' when destroying buffer... suggest re-writing app in Rust.`,
            );
        }
        this.pickbuffer = this.regl.framebuffer({
            colorFormat: "rgba",
            height: height,
            width: width,
        });

        //this is necessary, but I  don't know why?
        const loop = this.regl.frame(() => {});
        setTimeout(() => {
            this.refresh();
            loop.cancel();
        }, 200);
    }

    _getMousePosition(e) {
        const rect = this.canvas.getBoundingClientRect();
        return [e.clientX - rect.left, e.clientY - rect.top];
    }

    _getActualPosition(position) {
        const x = position[0] / this.x_scale - this.offset[0];
        const y = position[1] / this.y_scale - this.offset[1];
        return [x, y];
    }

    _getCanvasCoords(pos) {
        const x = (pos[0] + this.offset[0]) * this.x_scale;
        const y = (pos[1] + this.offset[1]) * this.y_scale;
        return [x, y];
    }

    getRange() {
        const tl = this._getActualPosition([0, 0]);
        const br = this._getActualPosition([this.width, this.height]);
        return {
            x_range: [tl[0], br[0]],
            y_range: [tl[1], br[1]],
            offset: [this.offset[0], this.offset[1]],
            scale: [this.x_scale, this.y_scale],
        };
    }

    setHighlightPoints(indexes) {
        if (indexes == null) {
            this.highlightPoints = null;
            return;
        }
        const hlp = {
            x_pos: [],
            y_pos: [],
            color: [],
            length: indexes.length,
            indexes: indexes,
        };
        if (this.mode === "3d") {
            hlp.z_pos = [];
        }
        for (const index of indexes) {
            hlp.x_pos.push(this.circles.x_pos[index]);
            hlp.y_pos.push(this.circles.y_pos[index]);
            if (this.mode === "3d") {
                hlp.z_pos.push(this.circles.z_pos[index]);
            }

            const st = index * 3;
            hlp.color.push(this.circles.color[st]);
            hlp.color.push(this.circles.color[st + 1]);
            hlp.color.push(this.circles.color[st + 2]);
        }
        this.highlightPoints = hlp;
    }

    setPointRadius(value) {
        if (!value || Number.isNaN(value)) {
            value = 0;
        }
        this.pointRadius = value;
    }
    setFilter(val) {
        this.isFiltered = val;
    }

    setHideOnFilter(value) {
        this.hideOnFilter = value;
    }

    setPointOpacity(value) {
        value = value > 1.0 ? 1.0 : value;
        value = value < 0.0 ? 0.0 : value;
        this.pointOpacity = value;
    }

    getPointOpacity() {
        return this.pointOpacity;
    }

    getPointRadius() {
        return this.pointRadius;
    }

    addLine(positionTo, positionFrom, color = [0, 0, 0], opacity = 1) {
        this.lines.position = this.lines.position.concat(
            positionTo,
            positionFrom,
        );
        this.lines.color = this.lines.color.concat(color, color);
        this.lines.opacity = this.lines.opacity.concat([opacity, opacity]);
        this.lines.count++;
        return this.lines.count - 1;
    }

    removeAllLines() {
        this.lines.position = [];
        this.lines.color = [];
        this.lines.opacity = [];
        this.lines.count = 0;
    }

    removeAllRectangles() {
        this.rects.position = [];
        this.rects.color = [];
        this.rects.opacity = [];
        this.rects.count = 0;
    }

    addRectangle(position, width, height, color = [0, 0, 0], opacity = 1) {
        this.rects.position.push(position);
        this.rects.position.push([position[0], position[1] + height]);
        this.rects.position.push([position[0] + width, position[1] + height]);
        this.rects.position.push([position[0] + width, position[1] + height]);
        this.rects.position.push([position[0] + width, position[1]]);
        this.rects.position.push([position[0], position[1]]);
        for (let a = 0; a < 6; a++) {
            this.rects.opacity.push(opacity);
        }
        const c = [color[0], color[1], color[2]];

        for (let a = 0; a < 6; a++) {
            this.rects.color.push(c);
        }
        this.rects.count++;
        return this.rects.count - 1;
    }

    changeImage(image, config, index) {
        this.images.props.text[index] = this.regl.texture({
            data: image,
            min: "linear",
        });
        this.resizeImage([config.width, config.height], index);
        this.setImagePosition(config.position[0], config.position[1], index);
    }

    changeImageOpacity(index, opacity) {
        const st = index * 6;
        const end = st + 6;
        for (let a = st; a < end; a++) {
            this.images.opacity[a] = opacity;
        }
    }

    moveImage(x, y) {
        const st = this.movingImage * 6;
        const en = st + 6;
        for (let n = st; n < en; n++) {
            const pos = this.images.position[n];
            pos[0] += x;
            pos[1] += y;
        }
        const xy = this.images.props.x_y[this.movingImage];
        xy[0] += x;
        xy[1] += y;
    }

    setImagePosition(x, y, index) {
        const p = this.images.position;
        const i = index * 6;
        const width = this.images.props.w_h[index][0];
        const height = this.images.props.w_h[index][1];
        y = -y - height;
        this.images.props.x_y[index] = [x, y];
        p[i] = [x, y];
        p[i + 1] = [x, y + height];
        p[i + 2] = [x + width, y + height];
        p[i + 3] = [x + width, y + height];
        p[i + 4] = [x + width, y];
        p[i + 5] = [x, y];
    }

    getImageDetails(index) {
        const wh = this.images.props.w_h[index];
        const xy = this.images.props.x_y[index];
        return {
            width: wh[0],
            height: wh[1],
            position: [xy[0], xy[1]],
            index: index,
        };
    }

    resizeImage(factor, index) {
        const st = index * 6;
        const wh = this.images.props.w_h[index];
        const xy = this.images.props.x_y[index];
        let width;
        let height = null;
        if (Array.isArray(factor)) {
            width = factor[0];
            height = factor[1];
            wh[0] = width;
            wh[1] = height;
        } else {
            width = wh[0] *= factor;
            height = wh[1] *= factor;
        }

        const x = xy[0];
        const y = xy[1];
        this.images.position[1] = [x, y + height];
        this.images.position[2] = [x + width, y + height];
        this.images.position[3] = [x + width, y + height];
        this.images.position[4] = [x + width, y];
    }

    addImage(image, config) {
        const c = config;

        const x = c.position[0];
        let y = -c.position[1];
        const height = c.height;
        const width = c.width;
        y = y - height;
        const image_index = this.images.position.length;
        this.images.position.push([x, y]);
        this.images.position.push([x, y + height]);
        this.images.position.push([x + width, y + height]);
        this.images.position.push([x + width, y + height]);
        this.images.position.push([x + width, y]);
        this.images.position.push([x, y]);
        this.image_position = [
            [-1, 1],
            [-1, 0],
            [0, 0],
            [0, 0],
            [0, 1],
            [-1.0],
        ];

        const opacity = c.opacity == null ? 1 : c.opacity;
        for (let a = 0; a < 6; a++) {
            this.images.color.push([1, 1, 1]);
            this.images.opacity.push(opacity);
        }
        this.images.count++;
        this.images.props.x_y.push([x, y]);
        this.images.props.w_h.push([width, height]);
        this.images.props.text.push(
            this.regl.texture({ data: image, min: "linear" }),
        );
        return this.images.count - 1;
    }

    removeImages() {
        this.images.position = [];
        this.images.color = [];
        this.images.opacity = [];
        this.images.props.x_y = [];
        this.images.props.w_h = [];
        this.images.props.text = [];
        this.images.count = 0;
    }

    addSquares(config) {
        this.squares.x_pos = config.x;
        this.squares.y_pos = config.y;
        const len = config.x.length;
        this.squares.localFilter = config.localFilter;
        this.squares.globalFilter = config.globalFilter;
        this.squares.pick_color = new Uint8Array(len * 3);
        this.squares.color = config.colors;
        for (let n = 0; n < len; n++) {
            const p = n * 3;
            const pb = this._getRGBFromIndex(n + 1);
            this.squares.pick_color[p] = pb[0];
            this.squares.pick_color[p + 1] = pb[1];
            this.squares.pick_color[p + 2] = pb[2];
        }
        this.squares.count = len;
    }

    /** this looks as though it will conveniently replace any existing circles, as opposed to adding
     * which is what I want... if not what the method name suggests.
     */
    addCircles(config) {
        this.circles.x_pos = config.x;
        this.circles.y_pos = config.y;
        if (this.mode === "3d") {
            this.circles.z_pos = config.z;
        }
        const len = config.x.length;
        this.circles.localFilter = config.localFilter;
        this.circles.globalFilter = config.globalFilter;
        this.circles.pick_color = new Uint8Array(len * 3);
        this.circles.color = new Uint8Array(len * 3);
        this.circles.colorFilter = new Uint8Array(len);
        for (let n = 0; n < len; n++) {
            const p = n * 3;
            const col = config.colorFunc(n);
            if (col) {
                this.circles.color[p] = col[0];
                this.circles.color[p + 1] = col[1];
                this.circles.color[p + 2] = col[2];
                this.circles.colorFilter[n] = 1;
            }
            const pb = this._getRGBFromIndex(n + 1) || [255, 0, 0];
            this.circles.pick_color[p] = pb[0];
            this.circles.pick_color[p + 1] = pb[1];
            this.circles.pick_color[p + 2] = pb[2];
        }
        this.circles.count = len;
    }

    updateSize(newSize, config) {
        this.circles.x_pos = config.x;
        this.circles.y_pos = config.y;
        if (this.mode === "3d") {
            this.circles.z_pos = config.z;
        }

        this.circles.localFilter = config.localFilter;
        this.circles.globalFilter = config.globalFilter;
        const newPickColor = new Uint8Array(newSize * 3);
        const newColor = new Uint8Array(newSize * 3);
        const newColorFilter = new Uint8Array(newSize);
        newColor.set(this.circles.color);
        newPickColor.set(this.circles.pick_color);
        newColorFilter.set(this.circles.colorFilter);
        this.circles.color = newColor;
        this.circles.pick_color = newPickColor;
        this.circles.colorFilter = newColorFilter;
        for (let n = this.circles.count; n < newSize; n++) {
            const p = n * 3;
            const col = config.colorFunc(n);
            if (col) {
                this.circles.color[p] = col[0];
                this.circles.color[p + 1] = col[1];
                this.circles.color[p + 2] = col[2];
                this.circles.colorFilter[n] = 1;
            }
            const pb = this._getRGBFromIndex(n + 1);
            this.circles.pick_color[p] = pb[0];
            this.circles.pick_color[p + 1] = pb[1];
            this.circles.pick_color[p + 2] = pb[2];
        }
        this.circles.count = newSize;
    }

    _getRGBFromIndex(index) {
        const b = Math.floor(index / 65536);
        const temp = index % 65536;
        const g = Math.floor(temp / 256);
        const r = temp % 256;
        return [r, g, b];
    }

    _getIndexFromRGB(rgb) {
        return rgb[2] * 65536 + rgb[1] * 256 + rgb[0];
    }

    _drawPickBuffer() {
        this.regl.clear({
            color: [0, 0, 0, 0],
            depth: 1,
            framebuffer: this.pickbuffer,
        });
        this._drawObjects(true);
    }
    //refesh all
    //in_view only those in view
    refresh() {
        //this.label_context.clearRect(0, 0, this.width, this.height);

        this.regl.clear({
            color: [0, 0, 0, 0],
            depth: 1,
        });
        this._drawObjects(false);
        this._drawPickBuffer();
        this.label_context.font = "30px Arial";
    }

    setCamera(distance, theta, phi) {
        if (this.mode !== "3d") {
            return;
        }
        const s = this.camera.cameraState;
        s.distance = Math.log(distance);
        s.theta = theta;
        s.phi = phi;
        this.camera.updateCamera();
    }

    getCameraSettings() {
        if (this.mode !== "3d") {
            return;
        }
        const s = this.camera.cameraState;
        return {
            distance: Math.E ** s.distance,
            theta: s.theta,
            phi: s.phi,
        };
    }

    zoom(amount) {
        this.x_scale *= amount;
        this.y_scale *= amount;
        this._drawObjects(false);
    }

    setFilterAction(type) {
        if (type === "grey") {
            this.hideOnFilter = false;
        } else {
            this.hideOnFilter = true;
        }
    }

    _drawObjects(buffer) {
        let obj = null;

        if (this.mode === "3d") {
            const cProj = this.camera.getProjection();
            obj = {
                cameraProjection: cProj.projection,
                cameraView: cProj.view,
                cameraDistance: cProj.distance,
                axisScales: this.axisScales,
                cameraEye: cProj.eye,
            };
        } else {
            obj = {
                x_scale: this.x_scale,
                y_scale: this.y_scale,
                offset: this.offset,
            };
        }
        obj.point_scale = this.pointScale
            ? (this.x_scale + this.y_scale) / 2 / this.pointScale
            : 1;
        obj.point_radius = this.pointRadius;
        obj.is_filtered = this.isFiltered ? 1 : 0;
        obj.point_opacity = this.pointOpacity;
        obj.hide_on_filter = this.hideOnFilter ? 1 : 0;

        for (const i of this.draw_order) {
            const type = this.object_types[i];
            //no objects of this type

            if (!type.data.count) {
                continue;
            }

            if (buffer) {
                if (!type.properties.pick_color) {
                    continue;
                }
                buffer = this.pickbuffer;
            } else {
                buffer = null;
            }

            obj.buffer = buffer;
            obj.count = type.data.count * type.vertices;
            obj.primitive = type.primitive;
            obj.is_buffer = buffer ? 1 : 0;

            //images special case - a draw commnad for each image
            if (i === 2) {
                for (let i = 0; i < type.data.count; i++) {
                    for (const prop in type.properties) {
                        obj[prop] = type.data[prop].slice(i * 6, i * 6 + 6);
                    }
                    for (const prop in this.images.props) {
                        obj[prop] = this.images.props[prop][i];
                    }
                    type.method(obj);
                }
                continue;
            }
            if (i === 3) {
                if (buffer) {
                    //continue;
                }

                if (this.x_scale >= this.y_scale) {
                    obj.side_length = 20 * this.x_scale;
                    obj.bottom_clip = this.y_scale / this.x_scale;
                    obj.right_clip = 1.0;
                } else {
                    obj.side_length = 20 * this.y_scale;
                    obj.right_clip = this.x_scale / this.y_scale;
                    obj.bottom_clip = 1.0;
                }
            }

            for (const prop in type.properties) {
                if (buffer) {
                    //swap color for pock buffer
                    if (prop === "pick_color") {
                        obj["color"] = type.data[prop];
                        continue;
                    }
                    if (prop === "color") {
                        continue;
                    }
                    obj[prop] = type.data[prop];
                } else {
                    if (prop === "pick_color") {
                        continue;
                    }
                    obj[prop] = type.data[prop];
                }
            }
            type.method(obj);
        }
        if (this.highlightPoints && !buffer) {
            obj.x_pos = this.highlightPoints.x_pos;
            obj.y_pos = this.highlightPoints.y_pos;
            obj.color = this.highlightPoints.color;
            obj.primitive = "points";
            obj.count = this.highlightPoints.length;
            if (this.mode === "3d") {
                obj.z_pos = this.highlightPoints.z_pos;
                this.__draw3DHighlights(obj);
            } else {
                this.__drawHighlights(obj);
            }
        }
    }

    getObjectsInRange(x, y, w, h) {
        console.log(`${x},${y},${w},${h}`);
        // there was a bug here when the x or y is -ve... patching it here
        if (x < 0) {
            w += x;
            x = 0;
        }
        if (y < 0) {
            h += y;
            y = 0;
        }
        // also an issue if the region extends beyond the canvas
        // pretty sure pickbuffer.width === this.width
        if (x + w > this.width) {
            w = this.width - x;
        }
        if (y + h > this.height) {
            h = this.height - y;
        }
        // Still have issues with dragging the region.
        // Could also happen elsewhere (ie with lasso tool? maybe not - does it use regl for picking?)
        // ideally, as a user, I'd like to be able to draw a region that goes off the canvas,
        // and have it still work - selecting objects that are off-screen.
        const max = w * h * 4;
        const pixels = this.regl.read({
            x: x,
            y: this.height - y - h,
            width: w,
            height: h,

            framebuffer: this.pickbuffer,
        });
        const s = new Set();
        for (let i = 0; i < max; i += 4) {
            const index = this._getIndexFromRGB([
                pixels[i],
                pixels[i + 1],
                pixels[i + 2],
            ]);
            if (index !== null) {
                s.add(index - 1);
            }
        }
        return s;
    }

    _getObjectAtPosition(position) {
        if (
            position[1] <= 0 ||
            position[0] <= 0 ||
            position[0] >= this.width ||
            position[1] >= this.height
        ) {
            return;
        }
        try {
            const pixel = this.regl.read({
                x: position[0],
                y: this.height - position[1],
                width: 1,
                height: 1,
                framebuffer: this.pickbuffer,
            });
            if (pixel[0] === 0 && pixel[1] === 0 && pixel[2] === 0) {
                return null;
            }
            //console.log(pixel);
            const index = this._getIndexFromRGB(pixel);
            if (index > 0) {
                return index - 1;
            }
            return null;
        } catch (e) {
            console.log(position);
            return null;
        }
    }

    addHandler(handler_type, handler, name) {
        const handler_dict = this.handlers[handler_type];
        if (!handler_dict) {
            throw "Handler Not Supported";
        }
        if (!name) {
            name = Object.keys(handler_dict).length;
        }
        handler_dict[name] = handler;
        return name;
    }

    removeHandler(handler_type, name) {
        const handler_dict = this.handlers[handler_type];
        if (!handler_dict) {
            throw "Handler Not Supported";
        }
        handler_dict["name"] = undefined;
    }

    _setUpBrush(origin) {
        const div = createEl(
            "div",
            {
                classes: ["wgl2di-brush"],
                styles: {
                    top: `${origin[1]}px`,
                    left: `${origin[0]}px`,
                },
            },
            this.div_container,
        );

        makeDraggable(div, {
            contain: true,
            ondragstart: () => (this.brush_moving = true),
            ondragend: () => this._brushingStopped(),
            doc: this.__doc__,
        });

        this.brush = { origin: origin, div: div, resizing: true };
    }

    _setUpPolyBrush(pos) {
        this.poly_brush = {
            points: [pos],
            active: true,
        };
        const ctx = this.label_context;
        ctx.strokeStyle = getComputedStyle(this.canvas).getPropertyValue(
            "color",
        );
        ctx.beginPath();
        ctx.moveTo(pos[0], pos[1]);
    }

    _extendPolyBrush(pos, threshold = 15) {
        const ctx = this.label_context;
        if (this.poly_brush.points.length) {
            const prev =
                this.poly_brush.points[this.poly_brush.points.length - 1];
            const dx = prev[0] - pos[0];
            const dy = prev[1] - pos[1];
            const d = Math.sqrt(dx ** 2 + dy ** 2);
            if (d < threshold) return;
        }

        ctx.lineTo(pos[0], pos[1]);
        ctx.stroke();
        this.poly_brush?.points.push(pos);
    }

    _finishPolyBrush() {
        if (this.poly_brush.points.length < 3) {
            this.clearBrush();
            return;
        }
        const ctx = this.label_context;
        const start = this.poly_brush.points[0];
        ctx.lineTo(start[0], start[1]);
        ctx.stroke();
        ctx.closePath();
        ctx.fillStyle = "lightgray";
        ctx.globalAlpha = 0.2;
        ctx.fill();
        ctx.globalAlpha = 1;
        const poly = [];
        this.poly_brush.active = false;
        for (const pt of this.poly_brush.points) {
            poly.push(this._getActualPosition(pt));
        }
        for (const i in this.handlers.brush_stopped) {
            setTimeout(() => this.handlers.brush_stopped[i](poly, true), 0);
        }
    }

    clearBrush() {
        if (this.brush) {
            this.brush.div.remove();
            this.brush = null;
        }
        if (this.poly_brush) {
            this.label_context.clearRect(0, 0, this.width, this.height);
            this.poly_brush = null;
        }
    }

    _brushingStopped() {
        if (!this.brush) return; //'cannot set properties of null' was possible here
        this.brush.resizing = false;
        const d = this.brush.div;
        const y = d.offsetTop;
        const x = d.offsetLeft;
        const w = d.offsetWidth;
        const h = d.offsetHeight;
        if (w < 2) {
            this.clearBrush();
            return;
        }
        if (this.mode === "3d") {
            const s = this.getObjectsInRange(x, y, w, h);

            for (const i in this.handlers.brush_stopped) {
                this.handlers.brush_stopped[i](s);
            }
        } else {
            const lt = this._getActualPosition([x, y]);
            const br = this._getActualPosition([x + w, y + h]);
            const info = {
                x_min: lt[0],
                x_max: br[0],
                y_max: -lt[1],
                y_min: -br[1],
            };
            for (const i in this.handlers.brush_stopped) {
                this.handlers.brush_stopped[i](info);
            }
        }
    }

    setGlobalFilter(filter) {
        this.circles.globalFilter = filter;
    }

    _finish(evt) {
        if (this.config.brush) {
            this.div_container.style.cursor = "crosshair";
        }
        if (this.brush?.resizing) {
            this._brushingStopped();
            //return;
        }
        if (this.poly_brush?.active) {
            this._finishPolyBrush();
        }

        //an object has finshed its drag
        if (this.object_clicked) {
            this.object_clicked = null;
            this.refresh(true);
        }

        //update which objects are now in view
        else {
            const ret = this.getRange();
            if (this.loop) {
                this.loop.cancel();
                this.loop = null;
            }

            //this.refresh();
            ret.imageMoved = this.imageMoved;
            this.imageMoved = null;
            for (const i in this.handlers.panning_stopped) {
                this.handlers.panning_stopped[i](ret);
            }
        }
        this.dragging = false;

        this.object_clicked = null;
        this.mouse_position = null;
    }

    _addHandlers() {
        const { noCameraControl } = this.config;
        //some listeners are on window while brush is active, and removed when brush is finished
        //allows dragging outside the canvas
        const mousemoveDragging = (e) => {
            if (this.brush) {
                if (this.brush.resizing) {
                    const origin = this.brush.origin;
                    const now = this._getMousePosition(e);
                    const left = `${Math.round(origin[0] < now[0] ? origin[0] : now[0])}px`;
                    const top = `${Math.round(origin[1] < now[1] ? origin[1] : now[1])}px`;
                    const width = `${Math.abs(origin[0] - now[0])}px`;
                    const height = `${Math.abs(origin[1] - now[1])}px`;
                    this.brush.div.style.top = top;
                    this.brush.div.style.left = left;
                    this.brush.div.style.width = width;
                    this.brush.div.style.height = height;

                    return;
                }
                if (this.brush.moving) {
                    this.dragging = false;
                    return;
                }
            }
            if (this.poly_brush?.active) {
                const pt = this._getMousePosition(e);
                this._extendPolyBrush(pt);
            }
            //is this a drag or just a click without the mouse moving
            if (this.mouse_position && !this.dragging) {
                const x_amount = e.pageX - this.mouse_position[0];
                const y_amount = e.pageY - this.mouse_position[1];
                if (Math.abs(x_amount) > 3 || Math.abs(y_amount) > 3) {
                    this.dragging = true;
                }
            }

            if (this.dragging) {
                const dx = e.pageX - this.mouse_position[0];
                const dy = e.pageY - this.mouse_position[1];

                if (this.mode === "3d") {
                    this.camera.mouseChange(dx / this.width, dy / this.height);
                } else {
                    const x_amount = dx / this.x_scale;
                    const y_amount = dy / this.y_scale;
                    if (
                        this.movingImage != null &&
                        e.ctrlKey &&
                        this.images.count > 0
                    ) {
                        this.moveImage(x_amount, y_amount);
                        this.imageMoved = this.getImageDetails(
                            this.movingImage,
                        );
                    } else if (this.offsets && e.which === 3) {
                        this.alterOffsets(x_amount, y_amount);
                    } else {
                        if (!this.config.lock_x_axis) {
                            this.offset[0] += x_amount;
                        }
                        this.offset[1] += y_amount;
                        if (!noCameraControl) {
                            for (const i in this.handlers.pan_or_zoom) {
                                this.handlers.pan_or_zoom[i](
                                    this.offset,
                                    this.x_scale,
                                    this.y_scale,
                                );
                            }
                        }
                    }
                }

                if (!this.loop) {
                    //self.label_context.clearRect(0, 0, self.width, self.height);
                    this.loop = this.regl.frame(() => {
                        this.regl.clear({
                            color: [0, 0, 0, 0],
                            depth: 1,
                        });
                        this._drawObjects(false);
                    });
                }
                this.mouse_position[1] = e.pageY;
                this.mouse_position[0] = e.pageX;
            }
        };
        this.div_container.addEventListener("mousemove", (e) => {
            if (this.dragging) return;
            //no drag event going on call any listners if mouse over/out an object
            const position = this._getMousePosition(e);
            const obj = this._getObjectAtPosition(position);
            if (obj != null && this.object_mouse_over == null) {
                for (const i in this.handlers["object_over"]) {
                    this.handlers.object_over[i](e, obj);
                }
                this.object_mouse_over = obj;
            } else if (obj == null && this.object_mouse_over) {
                for (const i in this.handlers["object_out"]) {
                    this.handlers.object_out[i](e, this.object_mouse_over);
                }
                this.object_mouse_over = null;
            }
            //move directly from one object to another
            else if (obj != null && obj !== this.object_mouse_over) {
                for (const i in this.handlers["object_over"]) {
                    this.handlers.object_over[i](e, obj);
                }
                this.object_mouse_over = obj;
            }
        });

        const mouseup = (evt) => {
            this.__doc__.removeEventListener("mousemove", mousemoveDragging);
            this.__doc__.removeEventListener("mouseup", mouseup);

            this._finish(evt);
            if (
                !this.dragging &&
                !(this.brush || this.poly_brush) &&
                evt.button === 0
            ) {
                const position = this._getMousePosition(evt);
                const obj = this._getObjectAtPosition(position);
                const apos = this._getActualPosition(position);
                if (obj) {
                    for (const i in this.handlers.object_clicked) {
                        this.handlers.object_clicked[i](obj, apos);
                    }
                }
            }
        };
        this.div_container.addEventListener(
            "wheel",
            (event) => {
                if (noCameraControl) return;
                event.preventDefault();
                const position = this._getActualPosition(
                    this._getMousePosition(event),
                );
                let wdelta = 0;
                if (event.wheelDelta > 0 || event.detail < 0) {
                    this.zoom_amount += 0.01;
                    wdelta = -30;
                } else {
                    this.zoom_amount -= 0.01;
                    wdelta = 30;
                }
                if (this.mode === "3d") {
                    this.camera.mouseWheel(0, wdelta);
                } else {
                    if (this.movingImage != null && event.ctrlKey) {
                        this.resizeImage(
                            1 + this.zoom_amount,
                            this.movingImage,
                        );
                        const new_position = this._getActualPosition(
                            this._getMousePosition(event),
                        );
                        this.moveImage(
                            new_position[0] - position[0],
                            new_position[1] - position[1],
                        );
                        this.imageMoved = this.getImageDetails(
                            this.movingImage,
                        );
                    } else {
                        if (!this.config.lock_x_axis) {
                            this.x_scale *= 1 + this.zoom_amount;
                        }
                        if (!this.config.lock_y_axis) {
                            this.y_scale *= 1 + this.zoom_amount;
                        }
                        const new_position = this._getActualPosition(
                            this._getMousePosition(event),
                        );
                        this.offset[0] += new_position[0] - position[0];
                        this.offset[1] += new_position[1] - position[1];
                        for (const i in this.handlers.pan_or_zoom) {
                            this.handlers.pan_or_zoom[i](
                                this.offset,
                                this.x_scale,
                                this.y_scale,
                            );
                        }
                    }
                }

                if (!this.loop) {
                    this.label_context.clearRect(0, 0, this.width, this.height);
                    this.loop = this.regl.frame(() => {
                        this.regl.clear({
                            color: [0, 0, 0, 0],
                            depth: 1,
                        });
                        this._drawObjects(false);
                    });
                }

                //clear the timeout user has not finished zooming
                clearTimeout(this._timer);
                //when user finishes call the esxpensive methods;
                this._timer = setTimeout(() => {
                    this.zoom_amount = 0;
                    this.loop.cancel();
                    this.loop = null;
                    this._drawPickBuffer(false);
                    if (this.brush) {
                    }
                    const ret = this.getRange();
                    ret.imageMoved = this.imageMoved;
                    this.imageMoved = null;
                    for (const name in this.handlers.zoom_stopped) {
                        this.handlers.zoom_stopped[name](ret);
                    }

                    this.refresh();
                }, 350);
            },
            { passive: false },
        );
        // https://github.com/Taylor-CCB-Group/MDV/issues/78 suggestion to add
        // {passive: true} above to help with Chrome issue on popout.
        // Doesn't fix it for me, although setting to *false* does get rid of a console error
        // "Unable to preventDefault inside passive event listener invocation."
        this.div_container.addEventListener("mousedown", (evt) => {
            this.__doc__.addEventListener("mousemove", mousemoveDragging);
            this.__doc__.addEventListener("mouseup", mouseup);

            if (evt.which === 3) {
                //add right click behaviour
            }
            const otherKey =
                evt.shiftKey ||
                evt.ctrlKey ||
                evt.which === 2 ||
                evt.which === 3;
            //create brush
            if (
                (this.config.brush && !otherKey) ||
                (!this.config.brush && otherKey)
            ) {
                const origin = this._getMousePosition(evt);
                if (this.config.brush === "default") {
                    const t = evt.target;
                    if (t.classList.contains("wgl2di-brush")) {
                        return;
                    }
                    if (this.brush) {
                        this.clearBrush();
                    }
                    const origin = this._getMousePosition(evt);
                    this._setUpBrush(origin);
                    return;
                }
                if (this.config.brush === "poly") {
                    if (this.poly_brush) {
                        this.clearBrush();
                    }
                    this._setUpPolyBrush(origin);
                    return;
                }
            }

            if (evt.otherKey) {
                this.div_container.style.cursor = "move";
            }

            const position = this._getMousePosition(evt);
            this.mouse_position = [evt.pageX, evt.pageY];
            this.clearBrush();
            evt.preventDefault();
        });
    }

    _initDrawMethods() {
        //loading images
        this.__drawCircles = this.regl({
            depth: this.draw_options.depth,
            blend: this.draw_options.blend,

            frag: `precision mediump float;
			varying vec3 fragColor;
			varying float op;
			varying float has_border;
			varying float grey_out;
			uniform float is_buffer;
			
			void main(){
			
				vec2 cxy = 2.0 * gl_PointCoord - 1.0;
                
				float r = dot(cxy, cxy);
				if (r > 1.0) {
					discard;
				}
				else{					
					if(r>0.60 && has_border==-1.0 && is_buffer==0.0 && grey_out==0.0){
					    gl_FragColor=vec4(0.1,0.1,0.1,1.0);
					}
					else{
						if (grey_out==1.0 && is_buffer==0.0){
							gl_FragColor = vec4(0.2,0.2,0.2,0.4);
						}
						else{
							gl_FragColor = vec4(fragColor,is_buffer==1.0?1.0:op);
						}
						
					}
				}
				
			}`,

            vert: `attribute float x_pos;
            attribute float y_pos;
			attribute vec3 color;
			attribute float lFilter;
			attribute float gFilter;
			attribute float colorFilter;
			varying vec3 fragColor;
			varying float op;
			varying float has_border;
			varying float grey_out;
			uniform float x_scale;
			uniform float y_scale;
			uniform vec2 offset;
			uniform float stage_height;
			uniform float stage_width;
			uniform float point_radius;
			uniform float is_filtered;
			uniform float point_opacity;
			uniform float hide_on_filter;
			uniform float point_scale;
			
			vec2 normalizeCoords(float posX, float posY){
				float x = (posX+offset[0])*x_scale;
				float y = (posY+offset[1])*y_scale;
				return vec2(2.0 * ((x / stage_width) - 0.5),-(2.0 * ((y / stage_height) - 0.5)));
			}
	
			void main() {
				float dd = 0.0;
				grey_out=0.0;
				if (colorFilter == 0.0) {
					return;
				}
				if (lFilter==2.0){
					return;
				}
			    if (gFilter>0.0 && gFilter != lFilter) {
					if(hide_on_filter==1.0){
						return;
					}
					else{
						grey_out=1.0;
					}
				}
				
			
				float r=point_radius;
			
				gl_PointSize = r*point_scale;
				vec3 c = vec3(255.0,255.0,255.0);
				if (color==c){
					op=0.1;
					fragColor=vec3(0.1,0.1,0.1);
				}
				else{
					op=point_opacity;
					fragColor = color/255.0;
				}
			
				
				
				has_border=lFilter-is_filtered;
			
				
				vec2 real_position = normalizeCoords(x_pos,-y_pos);
				gl_Position = vec4(real_position, 0.0, 1.0);
			}
			`,

            attributes: {
                x_pos: this.regl.prop("x_pos"),
                y_pos: this.regl.prop("y_pos"),
                color: this.regl.prop("color"),
                lFilter: this.regl.prop("localFilter"),
                gFilter: this.regl.prop("globalFilter"),
                colorFilter: this.regl.prop("colorFilter"),
            },

            uniforms: {
                x_scale: this.regl.prop("x_scale"),
                y_scale: this.regl.prop("y_scale"),
                stage_width: this.regl.context("viewportWidth"),
                stage_height: this.regl.context("viewportHeight"),
                offset: this.regl.prop("offset"),
                point_radius: this.regl.prop("point_radius"),
                point_opacity: this.regl.prop("point_opacity"),
                is_filtered: this.regl.prop("is_filtered"),
                is_buffer: this.regl.prop("is_buffer"),
                hide_on_filter: this.regl.prop("hide_on_filter"),
                point_scale: this.regl.prop("point_scale"),
            },

            count: this.regl.prop("count"),
            primitive: this.regl.prop("primitive"),
            framebuffer: this.regl.prop("buffer"),
        });

        this.__drawHighlights = this.regl({
            depth: this.draw_options.depth,
            blend: this.draw_options.blend,

            frag: `precision highp float;
			varying vec3 fragColor;
			
			void main(){

				vec2 cxy = 2.0 * gl_PointCoord - 1.0;
                
				float r = dot(cxy, cxy);
				if (r > 1.0) {
					discard;
				}
				else{					
					if(r>0.50) {
					    gl_FragColor=vec4(0.1,0.1,0.1,1.0);
					}
					else{
						gl_FragColor = vec4(fragColor,1);
					}
				}	
			}`,

            vert: `attribute float x_pos;
            attribute float y_pos;
			attribute vec3 color;
			uniform float point_radius;
			
			varying vec3 fragColor;
		
			uniform float x_scale;
			uniform float y_scale;
			uniform vec2 offset;
			uniform float stage_height;
			uniform float stage_width;
			uniform float point_scale;
			
			vec2 normalizeCoords(float posX, float posY){
				float x = (posX+offset[0])*x_scale;
				float y = (posY+offset[1])*y_scale;
				return vec2(2.0 * ((x / stage_width) - 0.5),-(2.0 * ((y / stage_height) - 0.5)));
			}
	
			void main() {
			
				float pr = point_radius*point_scale;
				pr = pr<10.0?10.0:pr;
				gl_PointSize = pr;
				fragColor = color/255.0;
					
				vec2 real_position = normalizeCoords(x_pos,-y_pos);
				gl_Position = vec4(real_position, 0.0, 1.0);
			}
			`,

            attributes: {
                x_pos: this.regl.prop("x_pos"),
                y_pos: this.regl.prop("y_pos"),
                color: this.regl.prop("color"),
            },

            uniforms: {
                x_scale: this.regl.prop("x_scale"),
                y_scale: this.regl.prop("y_scale"),
                stage_width: this.regl.context("viewportWidth"),
                stage_height: this.regl.context("viewportHeight"),
                offset: this.regl.prop("offset"),
                point_radius: this.regl.prop("point_radius"),
                point_scale: this.regl.prop("point_scale"),
            },

            count: this.regl.prop("count"),
            primitive: this.regl.prop("primitive"),
        });

        this.object_types[0]["method"] = this.__drawCircles;

        this.__drawLines = this.regl({
            // fragment shader
            frag: `precision highp float;
						varying vec3 fragColor;
						void main () {
							 gl_FragColor = vec4(fragColor,1);
						}`,

            vert: `
						attribute vec2 position;
						attribute vec3 color;
						attribute float opacity;
						uniform float x_scale;
						uniform float y_scale;
						uniform vec2 offset;
						uniform float stage_height;
						uniform float stage_width;
						varying vec3 fragColor;
						//varying float op;
						vec2 normalizeCoords(vec2 position){	    
							float x = (position[0]+offset[0])*x_scale;
							float y = (position[1]+offset[1])*y_scale;
				            return vec2(2.0 * ((x / stage_width) - 0.5),-(2.0 * ((y / stage_height) - 0.5)));
						}
						void main () {
							if (opacity==0.0){
								return;
							}
							fragColor=color/255.0;
							vec2 norm_pos =normalizeCoords(position);
							gl_Position = vec4(norm_pos, 0.0, 1.0);
						}`,
            attributes: {
                position: this.regl.prop("position"),
                color: this.regl.prop("color"),
                opacity: this.regl.prop("opacity"),
            },

            uniforms: {
                x_scale: this.regl.prop("x_scale"),
                y_scale: this.regl.prop("y_scale"),
                stage_height: this.regl.context("viewportHeight"),
                stage_width: this.regl.context("viewportWidth"),
                offset: this.regl.prop("offset"),
            },
            primitive: this.regl.prop("primitive"),
            count: this.regl.prop("count"),
        });
        this.object_types[1]["method"] = this.__drawLines;
        this.object_types[4]["method"] = this.__drawLines;

        this.__drawImages = this.regl({
            frag: `
				precision mediump float;
				uniform sampler2D text;
				varying vec2 uv;
				varying vec3 fragColor;
				varying float op;
				void main () {
					gl_FragColor = vec4(fragColor,op)*texture2D(text, uv);							
				}`,

            vert: `
				precision mediump float;

				attribute vec2 position;
				attribute vec3 color;
				attribute float opacity;

				uniform vec2 x_y;
				uniform vec2 w_h;
				uniform float stage_height;
				uniform float stage_width;
				uniform float x_scale;
				uniform float y_scale;
				uniform vec2 offset;

				varying vec2 uv;
				varying vec3 fragColor;
				varying float op;
		
				
				vec2 normalizeCoords(vec2 pos){
					float x = (pos[0]+offset[0])*x_scale;
					float y = (pos[1]+offset[1])*y_scale;
					return vec2(2.0 * ((x / stage_width) - 0.5),-(2.0 * ((y / stage_height) - 0.5)));
				}

				void main () {
					if (opacity==0.0){
						return;
					}
					op=opacity;
					
					vec2 new_pos=normalizeCoords(position);
					fragColor = color;
				
					float x_factor = 1.0/(((w_h[0]*x_scale)/stage_width)*2.0);
					float y_factor = 1.0/(((w_h[1]*y_scale)/stage_height)*2.0);
					

					float x_offset=(((x_y[0]+offset[0])*x_scale)/stage_width)*2.0*x_factor;
					float y_offset=(((x_y[1]+offset[1])*y_scale)/stage_height)*2.0*y_factor;
					uv = vec2((new_pos[0]*x_factor)+x_factor-x_offset,-(new_pos[1]*y_factor)+y_factor-y_offset);
					gl_Position = vec4(new_pos, 0.0, 1.0);

				}`,

            attributes: {
                position: this.regl.prop("position"),
                color: this.regl.prop("color"),
                opacity: this.regl.prop("opacity"),
            },

            uniforms: {
                stage_height: this.regl.context("viewportHeight"),
                stage_width: this.regl.context("viewportWidth"),
                w_h: this.regl.prop("w_h"),
                x_y: this.regl.prop("x_y"),
                text: this.regl.prop("text"),

                x_scale: this.regl.prop("x_scale"),
                y_scale: this.regl.prop("y_scale"),
                offset: this.regl.prop("offset"),
            },

            count: this.regl.prop("count"),
            depth: this.draw_options.depth,
            blend: this.draw_options.blend,
        });

        this.object_types[2]["method"] = this.__drawImages;

        this.__drawSquares = this.regl({
            frag: `
				precision highp float;
				varying vec3 fragColor;
				varying float r_clip;
				varying float b_clip;
				uniform int is_buffer;
				void main(){
					if (gl_PointCoord[1]<1.0-b_clip){
						discard;
					}
					if (gl_PointCoord[0]> r_clip){
						discard;
					}				
					gl_FragColor = vec4(fragColor,1);					
				}
				`,

            vert: ` 
				attribute float x_pos;
				attribute float y_pos;			
				attribute vec3 color;
			
				varying float r_clip;
				varying float b_clip;
				varying vec3 fragColor;

				uniform float x_scale;
				uniform float y_scale;
				uniform vec2 offset;
				uniform float stage_height;
				uniform float stage_width;
				uniform float right_clip;
				uniform float bottom_clip;
				uniform float side_length;
				
				vec2 normalizeCoords(float posX,float posY){
					float x = (posX+offset[0])*x_scale;
					float y = (posY+offset[1])*y_scale;
					return vec2(2.0 * ((x / stage_width) - 0.5),-(2.0 * ((y / stage_height) - 0.5)));
				}
				void main() {
					gl_PointSize = side_length;
					fragColor = color/255.0;
					r_clip=right_clip;
					b_clip= bottom_clip;
					float y = y_pos;
					float x = x_pos;
					if (bottom_clip!=1.0){
						y-=10.0*(1.0/b_clip)-10.0;
					}
					if (r_clip!=1.0){
						x+=10.0*(1.0/r_clip)-10.0;
					}

					vec2 real_position = normalizeCoords(x,y);
				
					gl_Position = vec4(real_position, 0.0, 1.0);
				}
			`,

            attributes: {
                x_pos: this.regl.prop("x_pos"),
                y_pos: this.regl.prop("y_pos"),
                color: this.regl.prop("color"),
            },

            uniforms: {
                x_scale: this.regl.prop("x_scale"),
                y_scale: this.regl.prop("y_scale"),
                stage_height: this.regl.context("viewportHeight"),
                stage_width: this.regl.context("viewportWidth"),
                offset: this.regl.prop("offset"),
                is_buffer: this.regl.prop("is_buffer"),
                right_clip: this.regl.prop("right_clip"),
                bottom_clip: this.regl.prop("bottom_clip"),
                side_length: this.regl.prop("side_length"),
            },

            count: this.regl.prop("count"),
            primitive: this.regl.prop("primitive"),
            framebuffer: this.regl.prop("buffer"),
        });
        this.object_types[3]["method"] = this.__drawSquares;

        this.__draw3DCircles = this.regl({
            frag: `
            precision highp float;
            varying vec3 fragColor;
			varying float op;
			varying float has_border;
			uniform float is_buffer;
			varying float grey_out;
            void main () {
			  float r = length(gl_PointCoord.xy - 0.4) ;
              if (r> 0.4) {
                discard;
              }
			  else{
				if(r>0.3 && has_border==-1.0 && is_buffer==0.0 && grey_out==0.0){
					gl_FragColor=vec4(0.1,0.1,0.1,1.0);
				}
				else{
					if (grey_out==1.0 && is_buffer==0.0){
						gl_FragColor = vec4(0.2,0.2,0.2,0.4);
					}
					else{
						gl_FragColor = vec4(fragColor,is_buffer==1.0?1.0:op);
					}
					
				}

			  }
             
            }
            `,

            vert: `
            precision highp float;

            attribute float x_pos;
            attribute float y_pos;
            attribute float z_pos;
            attribute vec3 color;
            attribute float lFilter;
            attribute float gFilter;
            attribute float colorFilter;
			
			uniform float time;
            uniform float cdistance;
            uniform mat4 view, projection;
			uniform vec3 eye;
			uniform float point_radius;
			uniform float point_opacity;
			uniform float is_buffer;
			uniform float is_filtered;
			uniform float hide_on_filter;
			uniform vec3 axis_scales;


            varying vec3 fragColor;
			varying float op;
			varying float has_border;
			varying float grey_out;

            void main () {
				grey_out=0.0;
				if (colorFilter == 0.0) {
					return;
				}
			    if (gFilter>0.0 && gFilter != lFilter) {
					if(hide_on_filter==1.0){
						return;
					}
					else{
						grey_out=1.0;
					}
				}
                vec3 position = vec3(-x_pos *axis_scales.x,y_pos*axis_scales.y,-z_pos*axis_scales.z);
                fragColor=color/255.0;
				op = is_buffer==1.0?1.0:point_opacity;
				has_border=lFilter-is_filtered;
               
				gl_PointSize = point_radius; //(distance(eye, position.xyz) /30.0);
                gl_Position = projection * view * vec4(2.0 * position, 1.0);
            }
            `,

            attributes: {
                x_pos: this.regl.prop("x_pos"),
                y_pos: this.regl.prop("y_pos"),
                z_pos: this.regl.prop("z_pos"),
                color: this.regl.prop("color"),
                lFilter: this.regl.prop("localFilter"),
                gFilter: this.regl.prop("globalFilter"),
                colorFilter: this.regl.prop("colorFilter"),
            },

            uniforms: {
                point_radius: this.regl.prop("point_radius"),
                point_opacity: this.regl.prop("point_opacity"),
                is_buffer: this.regl.prop("is_buffer"),
                view: this.regl.prop("cameraView"),
                cdistance: this.regl.prop("cameraDistance"),
                eye: this.regl.prop("cameraEye"),
                is_filtered: this.regl.prop("is_filtered"),
                hide_on_filter: this.regl.prop("hide_on_filter"),
                axis_scales: this.regl.prop("axisScales"),

                projection: (context, props) => {
                    return glMatrix.mat4.perspective(
                        props.cameraProjection,
                        Math.PI / 4,
                        context.viewportWidth / context.viewportHeight,
                        1,
                        100000,
                    );
                },
            },

            count: this.regl.prop("count"),
            framebuffer: this.regl.prop("buffer"),

            blend: this.draw_options.blend,
            primitive: "points",
        });

        this.__draw3DHighlights = this.regl({
            frag: `
            precision mediump float;
            varying vec3 fragColor;
		
            void main () {
			  float r = length(gl_PointCoord.xy - 0.4) ;
              if (r> 0.4) {
                discard;
              }
			  else{
				if(r>0.3){
					gl_FragColor=vec4(0.1,0.1,0.1,1.0);
				}
				else{
					
						gl_FragColor = vec4(fragColor,1.0);
				}
					
			  }
             
            }
            `,

            vert: `
            precision mediump float;

            attribute float x_pos;
            attribute float y_pos;
            attribute float z_pos;
            attribute vec3 color;
			
			uniform float time;
    
            uniform mat4 view, projection;
			uniform vec3 eye;
			uniform float point_radius;
			uniform vec3 axis_scales;
		


            varying vec3 fragColor;
		

            void main () {
			
                vec3 position = vec3(-x_pos*axis_scales.x,y_pos*axis_scales.y,-z_pos*axis_scales.z);
                fragColor=color/255.0;
               
				gl_PointSize = point_radius*2.0; //(distance(eye, position.xyz) /30.0);
                gl_Position = projection * view * vec4(2.0 * position, 1.0);
            }
            `,

            attributes: {
                x_pos: this.regl.prop("x_pos"),
                y_pos: this.regl.prop("y_pos"),
                z_pos: this.regl.prop("z_pos"),
                color: this.regl.prop("color"),
            },

            uniforms: {
                point_radius: this.regl.prop("point_radius"),

                view: this.regl.prop("cameraView"),
                cdistance: this.regl.prop("cameraDistance"),
                eye: this.regl.prop("cameraEye"),
                axis_scales: this.regl.prop("axisScales"),

                projection: (context, props) => {
                    return glMatrix.mat4.perspective(
                        props.cameraProjection,
                        Math.PI / 4,
                        context.viewportWidth / context.viewportHeight,
                        1,
                        100000,
                    );
                },
            },

            count: this.regl.prop("count"),
            blend: this.draw_options.blend,
            primitive: "points",
        });

        this.__draw3DLines = this.regl({
            // fragment shader
            frag: " precision highp float;\n\
                    varying vec3 fragColor;\n\
                    void main () {\n\
                         gl_FragColor = vec4(fragColor,1);\n\
                    }\n",

            vert: `
                    attribute vec3 position;
                    attribute vec3 color;
                    uniform mat4 view, projection;
					uniform vec3 axis_scales;
    
                    varying vec3 fragColor;
                  
                    
                    void main () {
                     
                        fragColor=color;
                        gl_Position = projection * view * vec4(2.0 * position*axis_scales*vec3(-1.0,1.0,-1.0), 1.0);
                       
                       
                    }`,
            attributes: {
                position: this.regl.prop("position"),
                color: this.regl.prop("color"),
            },
            uniforms: {
                view: this.regl.prop("cameraView"),
                axis_scales: this.regl.prop("axisScales"),
                projection: (context, props) => {
                    return glMatrix.mat4.perspective(
                        props.cameraProjection,
                        Math.PI / 4,
                        context.viewportWidth / context.viewportHeight,
                        0.01,
                        100000,
                    );
                },
            },

            primitive: this.regl.prop("primitive"),
            framebuffer: this.regl.prop("buffer"),
            count: this.regl.prop("count"),
            depth: { enable: false },
        });
    }
}

export { WGL2DI };
