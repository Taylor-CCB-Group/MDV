import { getProjectURL } from "../dataloaders/DataLoaderUtil.ts";
import { createEl } from "../utilities/Elements.js";

class ImageTable {
    /**
     * @param {HTMLElement|string} parent_div id or element
     * @param {DataStore} data_view
     * @param {Object} config
     */
    constructor(parent_div, data_view, config) {
        if (typeof parent_div === "string") {
            parent_div = document.getElementById(`#${parent_div}`);
        }
        if (!config) {
            config = {};
        }
        this.domId = this._getRandomString();
        this.row_first = 0;
        this.row_last = 139;
        this.tile_height = 205;
        this.tile_width = 205;
        this.__doc__ = document;
        this.display_columns = [];
        this.selection_mode = true;
        this.imageFunc = config.imageFunc;
        //The base url contains the stub of the image e.g. /images/img
        //and is not just the folder- need to remove trailing slash
        //Although clunky this allows integer columns to represent the img
        //eg. 1 becomes /images/im1.png
        
        this.base_url = getProjectURL(config.base_url)
        if (!config.base_url.endsWith("/")) {
            this.base_url = this.base_url.replace(/\/+$/, "");
        }
        
      
        this.parent = parent_div;
        this.selected_tiles = {};

        this.cache_size = 5;
        const bb = this.parent.getBoundingClientRect();
        this.config = config;
        const bg = config.background_color;
        this.view_port = createEl(
            "div",
            {
                styles: {
                    height: `${bb.height}px`,
                    width: `${bb.width}px`,
                    overflow: "auto",
                    // backgroundColor:bg
                },
            },
            this.parent,
        );
        this.canvas = createEl(
            "div",
            {
                styles: {
                    position: "relative",
                },
            },
            this.view_port,
        );
        this.canvas.addEventListener("click", (e) => {
            this.imageClicked(e, e.srcElement);
        });

        this.canvas.addEventListener("mouseover", (e) => {
            this.mouseOver(e, e.srcElement);
        });

        this.canvas.addEventListener("keydown", (e) => {
            this.keyPressed(e);
        });

        this.data_view = data_view;
        this._setCanvasHeight();
        this.view_port.addEventListener("scroll", (e) => this._hasScrolled(e));
        this.parent.append(this.view_port);
        /*let end_row = Math.floor((this.height+(this.cache_size*this.tile_height))/this.tile_height);
        this.max_difference= end_row+this.cache_size;*/
        // this._addListeners();
        //this.render(0,end_row,true);
        this.resize_timeout = null;
        this.resize_timeout_length = 50;
        this.listeners = {
            image_clicked: new Map(),
            image_over: new Map(),
            image_out: new Map(),
            data_changed: new Map(),
            image_selected: new Map(),
        };
        this.highlightcolors = null;
        this.tag_color_palletes = {};
        this.image_suffix = config.image_suffix ? config.image_suffix : ".png";

        if (config.columns) {
            this.setColumns(config.columns);
        }
        //load first image to get the orignal dimensions
        const im = new Image();
        im.crossOrigin = "anonymous";
        im.onload = (e) => {
            this.originalDimensions = [im.width, im.height];
            this.preferred_width = this.img_width = im.width;
            this.preferred_height = this.img_height = im.height;
            if (config.initial_image_width) {
                this.preferred_width = config.initial_image_width;
                this.preferred_height = Math.round(
                    (this.preferred_width / this.img_width) * this.img_height,
                );
            }
            this._resize();
        };
        //load in first image to get preferred dimensions
        //if none are given
        const type = this.config.image_type;
        const image = this.data_view.getItemField(1, this.config.image_key);
        im.src = `${this.base_url}${image}.${type}`;
    }

    mouseOver(e, img) {
        if (img.id) {
            const arr = img.id.split("-");
            const index = Number.parseInt(arr[1]);
            const item = this.data_view.getItem(index);
            if (!item) {
                return;
            }
            if (this.image_over) {
                //already in image
                if (item.id === this.image_over.id) {
                    return;
                }
                //move from one image to another
                this.listeners.image_out.forEach((func) => {
                    func(e, this.image_over[0], self.image_over[1]);
                });
            }
            this.listeners.image_over.forEach((func) => {
                func(e, item, img);
            });
            this.image_over = img;
        }
        //mouse out
        else {
            if (this.image_over) {
                this.listeners.image_out.forEach((func) => {
                    func(e, self.image_over[0], self.image_over[1]);
                });
            }
            this.image_over = null;
        }
    }

    imageClicked(e, img) {
        let id = img.id;
        if (!id) {
            id = img.parentElement.id;
        }
        if (!id) {
            return;
        }
        const arr = id.split("-");
        if (arr[0] === "mlvtile") {
            let ids = null;

            const index = Number.parseInt(arr[2]);
            if (e.shiftKey && self.last_img_clicked != null) {
                const l_index = this.last_img_clicked.index;

                range = [];
                ids = [];
                const diff = index - l_index < 0 ? -1 : 1;
                let st = l_index + 1;
                let en = index + 1;
                if (diff === -1) {
                    st = index;
                    en = l_index;
                }
                for (let i = st; i < en; i++) {
                    const id = this.data_view.getId(i);

                    ids.push(id);
                }
                ids.push(this.data_view.getId(l_index));
            }

            const id = this.data_view.getId(index);

            this.listeners.image_clicked.forEach((func) => {
                func(e, id, img);
            });

            this.last_index_clicked = index;

            if (this.selection_mode) {
                this.setSelectedTiles([index], e.ctrlKey, true);
            }
        }
    }

    keyPressed(e) {
        if (e.which === 39 || e.which === 40) {
            if (
                this.last_index_clicked ===
                this.data_view.getFilteredItems().length - 1
            ) {
                return;
            }
            if (this.last_index_clicked || this.last_index_clicked === 0) {
                this.last_index_clicked++;
                this.setSelectedTiles(
                    [this.data_view.getItem(this.last_index_clicked).id],
                    false,
                    true,
                );
                if (this.last_index_clicked > this.getLastTileInView()) {
                    this.show(this.last_index_clicked);
                }
            }
        } else if (e.which === 37 || e.which === 38) {
            if (this.last_index_clicked === 0) {
                return;
            }
            if (this.last_index_clicked) {
                this.last_index_clicked--;
                this.setSelectedTiles(
                    [this.data_view.getItem(this.last_index_clicked).id],
                    false,
                    true,
                );
                if (this.last_index_clicked < this.getFirstTileInView()) {
                    this.show(this.last_index_clicked);
                }
            }
        }
    }

    addListener(type, func, id) {
        const listener = this.listeners[type];
        if (!listener) {
            return null;
        }
        if (!id) {
            id = this._getRandomString();
        }
        listener.set(id, func);
        return id;
    }

    removeListener(type, id) {
        const listener = this.listeners[type];
        if (!listener) {
            return false;
        }
        return listener.delete(id);
    }

    setImageWidth(width, redraw) {
        this.setImageDimensions(
            [width, Math.round((width / this.img_width) * this.img_height)],
            redraw,
        );
    }

    setImageLabel(col) {
        this.config.image_label = col;
        const fti = this.getFirstTileInView();
        this.setImageDimensions();
        this.show(fti);
    }
    setImageTitle(col) {
        this.config.image_title = col;
        const fti = this.getFirstTileInView();
        this.show(fti);
    }

    setImageDimensions(dim, redraw) {
        this.lrm = this.config.margins.left_right;
        this.tbm = this.config.margins.top_bottom;
        const ft = this.getFirstTileInView();
        if (!dim) {
            dim = [this.preferred_width, this.preferred_height];
        } else {
            this.preferred_width = dim[0];
            this.preferred_height = dim[1];
        }

        if (
            this.width < this.preferred_width + this.lrm * 2 &&
            this.width !== 0
        ) {
            dim[0] = this.width - this.lrm * 3;
            this.num_per_row = 1;
        } else {
            this.num_per_row = Math.ceil(
                (this.width - this.lrm) / (this.preferred_width + this.lrm),
            );
            dim[0] =
                Math.floor((this.width - 2 * this.lrm) / this.num_per_row) -
                this.lrm;
        }
        dim[1] = Math.round(
            (dim[0] / this.preferred_width) * this.preferred_height,
        );

        if (this.config.image_label) {
            this.label_size = Math.round(dim[1] / 12);
            this.label_size = Math.max(this.label_size, 12);
            //need to fit label in
            this.tbm =
                this.tbm < this.label_size + 4 ? this.label_size + 4 : this.tbm;
        }

        this.tile_width = Number.parseInt(dim[0]) + this.lrm;
        this.tile_height = Number.parseInt(dim[1]) + this.tbm;
        this.t_width = Number.parseInt(dim[0]);
        this.t_height = Number.parseInt(dim[1]);
        const end_row = Math.floor(
            (this.height + this.cache_size * this.tile_height) /
                this.tile_height,
        );
        this.max_difference = end_row + this.cache_size;
        if (redraw) {
            this.show(ft);
        }
    }
    setPixelated(pixelated) {
        this.view_port.style.imageRendering = pixelated ? "pixelated" : "auto";
    }

    resize() {
        clearTimeout(this.resize_timeout);
        this.resize_timeout = setTimeout(() => {
            this._resize();
        }, this.resize_timeout_length);
    }

    getFirstTileInView() {
        const top = this.view_port.scrollTop;
        return Math.floor(top / this.tile_height) * this.num_per_row;
    }

    getLastTileInView() {
        const bottom = this.view_port.scrollTop + this.view_port.height();
        return Math.floor(bottom / this.tile_height) * this.num_per_row;
    }

    setColumns(columns) {
        this.columns = columns;
    }

    setColorBy(color_by, overlay) {
        this.color_by = color_by;
        this.color_overlay = overlay;
    }

    _setCanvasHeight() {
        const h =
            Math.ceil(this.data_view.getLength() / this.num_per_row) *
                this.tile_height +
            this.tbm;
        this.canvas.style.height = `${h}px`;
    }

    _hasScrolled() {
        clearTimeout(this.scroll_timeout);
        const s_top = this.view_port.scrollTop;
        let begin_row = Math.floor(
            (s_top - this.cache_size * this.tile_height) / this.tile_height,
        );
        const end_row = Math.floor(
            (s_top + this.height + this.cache_size * this.tile_height) /
                this.tile_height,
        );
        if (begin_row < 0) {
            begin_row = 0;
        }
        let elapse = 10;
        if (
            Math.abs(begin_row - this.row_displayed_first) > this.max_difference
        ) {
            elapse = 50;
        }
        this.scroll_timeout = setTimeout(() => {
            if (
                Math.abs(begin_row - this.row_displayed_first) >
                this.max_difference
            ) {
                this.render(begin_row, end_row, true);
            } else {
                this.render(begin_row, end_row);
            }
        }, elapse);
    }

    _removeByClass(className) {
        document.querySelectorAll(className).forEach((e) => e.remove());
    }

    render(begin_row, end_row, all) {
        if (all) {
            this.canvas.innerHTML = "";
            for (let n = begin_row; n < end_row; n++) {
                this._addRow(n);
            }
        } else if (begin_row === this.row_displayed_first && begin_row !== 0) {
            return;
        } else {
            if (begin_row < this.row_displayed_first) {
                for (let n = begin_row; n < this.row_displayed_first; n++) {
                    this._addRow(n);
                }
                for (let n = end_row; n < this.row_displayed_last; n++) {
                    this._removeByClass(`.mlv-tile-${this.domId}-${n}`);
                }
            } else {
                for (let n = this.row_displayed_last; n < end_row; n++) {
                    this._addRow(n);
                }

                for (let n = this.row_displayed_first; n < begin_row; n++) {
                    this._removeByClass(`.mlv-tile-${this.domId}-${n}`);
                }
            }
            if (begin_row === 0) {
                for (let n = end_row; n < this.row_displayed_last; n++) {
                    this._removeByClass(`.mlv-tile-${this.domId}-${n}`);
                }
            }
        }
        this.row_displayed_first = begin_row;
        this.row_displayed_last = end_row;
    }

    _resize(fti) {
        if (!fti) {
            fti = this.getFirstTileInView();
        }
        const bb = this.parent.getBoundingClientRect();
        this.width = Math.round(bb.width);
        this.height = Math.round(bb.height);
        this.view_port.style.height = `${this.height}px`;
        this.view_port.style.width = `${this.width}px`;
        this.setImageDimensions();
        this.show(fti);
    }

    _calculateTopBottomRow(first_tile_index) {
        if (!first_tile_index) {
            first_tile_index = 0;
        }
        let s_top = 0;
        if (first_tile_index || first_tile_index === 0) {
            s_top =
                Math.floor(first_tile_index / this.num_per_row) *
                this.tile_height;
        } else {
            s_top = this.view_port.scrollTop;
        }
        const bb = this.view_port.getBoundingClientRect();
        const height = bb.height;
        let begin_row = Math.floor(
            (s_top - this.cache_size * this.tile_height) / this.tile_height,
        );
        const end_row = Math.floor(
            (s_top + height + this.cache_size * this.tile_height) /
                this.tile_height,
        );
        if (begin_row < 0) {
            begin_row = 0;
        }
        /*if (begin_row==0){
            s_top=0;
        }*/
        return { top: begin_row, bottom: end_row, scroll_top: s_top };
    }

    setSelectedTiles(ids, append, propagate) {
        if (!append) {
            for (const id in this.selected_tiles) {
                const t = this.__doc__.getElementById(
                    `mlvtile-${this.domId}-${id}`,
                );
                if (t) {
                    t.style.border = this.selected_tiles[id];
                }
            }
            this.selected_tiles = {};
        }
        for (const id of ids) {
            const tile = this.__doc__.getElementById(
                `mlvtile-${this.domId}-${id}`,
            );
            if (tile) {
                this.selected_tiles[id] = tile.style.border;
                tile.style.border = "4px solid goldenrod";
                console.log(tile.style.border);
            }
        }
        if (propagate) {
            this.listeners.image_selected.forEach((func) => {
                func(ids);
            });
        }
    }

    scrollToTile(index, select) {
        const pos = this.data_view.data.indexOf(index);
        const obj = this._calculateTopBottomRow(pos);
        if (
            Math.abs(obj.top - this.row_displayed_first) > this.max_difference
        ) {
            this.render(obj.top, obj.bottom, true);
        } else {
            this.render(obj.top, obj.bottom);
        }
        this.view_port.scrollTop = obj.scroll_top;
        if (select) {
            this.setSelectedTiles([pos]);
        }
    }

    _getRandomString(len, an) {
        if (!len) {
            len = 6;
        }
        an = an?.toLowerCase();
        let str = "";
        let i = 0;
        const min = an === "a" ? 10 : 0;
        const max = an === "n" ? 10 : 62;
        while (i++ < len) {
            let r = (Math.random() * (max - min) + min) << 0;
            str += String.fromCharCode((r += r > 9 ? (r < 36 ? 55 : 61) : 48));
        }
        return str;
    }

    show(first_tile_index) {
        if (!first_tile_index) {
            first_tile_index = this.getFirstTileInView();
        }
        this._setCanvasHeight();
        const obj = this._calculateTopBottomRow(first_tile_index);
        this.render(obj.top, obj.bottom, true);
        this.view_port.scrollTop = obj.scroll_top;
    }

    _addRow(row) {
        const label = this.config.image_label;
        const st = row * this.num_per_row;
        const en = st + this.num_per_row;
        const top = row * this.tile_height + this.tbm;
        let x = 0;
        const w = `${this.t_width}px`;
        const h = `${this.t_height}px`;
        const type = this.config.image_type;
        // looking for something to provide a default name; "Gene" is ok for ytrap...
        const titleColumn = this.config.image_title || "Gene";
        const hasTitleColumn =
            this.data_view.dataStore.columnIndex[titleColumn] !== undefined;
        for (let i = st; i < en; i++) {
            const id = this.data_view.getId(i);
            if (id == null) {
                return;
            }
            const image = this.data_view.getItemField(i, this.config.image_key);
            const text = hasTitleColumn
                ? `${titleColumn} '${this.data_view.getItemField(i, titleColumn)}'`
                : image;
            const missing = image === "missing";
            const left = x * this.tile_width + this.lrm;
            x++;
            let border = "";
            if (this.color_by) {
                const color = this.color_by(id);
                border = `4px solid ${color}`;
                if (this.color_overlay) {
                    createEl(
                        "div",
                        {
                            styles: {
                                height: h,
                                width: w,
                                left: `${left}px`,
                                top: `${top}px`,
                                zIndex: 1,
                                pointerEvents: "none",
                                backgroundColor: color,
                                opacity: this.color_overlay,
                            },
                            classes: [
                                "mlv-tile",
                                `mlv-tile-overlay-${this.domId}`,
                                `mlv-tile-${this.domId}-${row}`,
                            ], //this class used to remove rows no longer in view
                        },
                        this.canvas,
                    );
                }
            }
            const extra_classes = [];
            if (this.selected_tiles[id]) {
                border = "4px solid goldenrod";
            }
            if (missing) {
                extra_classes.push("mlv-tile-missing");
            }
            createEl(
                "img",
                {
                    src: `${this.base_url}${image}.${type}`,
                    id: `mlvtile-${this.domId}-${i}`,
                    styles: {
                        border: border,
                        height: h,
                        width: w,
                        left: `${left}px`,
                        top: `${top}px`,
                    },
                    alt: text,
                    title: text,
                    crossorigin: "",
                    classes: [
                        "mlv-tile",
                        `mlv-tile-${this.domId}`,
                        `mlv-tile-${this.domId}-${row}`,
                    ].concat(extra_classes),
                },
                this.canvas,
            );
            if (label) {
                createEl(
                    "div",
                    {
                        text: this.data_view.getItemField(i, label),
                        classes: ["mdv-image-label"],
                        styles: {
                            width: w,
                            fontSize: `${this.label_size}px`,
                            top: `${top - this.label_size}px`,
                            left: `${left}px`,
                        },
                    },
                    this.canvas,
                );
            }
        }
    }
}

export default ImageTable;
