import Split from "split.js";
import { createEl } from "./ElementsTyped";

//PJT: moved to ElementsTyped.ts because I realised type annotations broke jsdoc script
//(even though there were understood correctly by language server / IDE)
//function createEl(type,attrs,parent){

function createSVGEl(type, attrs, parent) {
    const el = document.createElementNS("http://www.w3.org/2000/svg", type);

    if (attrs) {
        for (const idx in attrs) {
            if ((idx === "styles" || idx === "style") && typeof attrs[idx] === "object") {
                for (const prop in attrs[idx]) {
                    el.style[prop] = attrs[idx][prop];
                }
            } else if (idx === "text") {
                el.textContent = attrs[idx];
            } else if (idx === "classes") {
                for (const cl of attrs[idx]) {
                    el.classList.add(cl);
                }
            } else {
                el.setAttributeNS(null, idx, attrs[idx]);
            }
        }
    }
    if (parent) {
        parent.append(el);
    }
    return el;
}

function splitPane(el, config = {}) {
    const dir = config.direction || "horizontal";
    const number = config.number || 2;
    const panes = [];
    const classes = dir === "horizontal" ? ["split-horizontal"] : ["split-vertical"];
    for (let i = 0; i < number; i++) {
        panes.push(createEl("div", { classes: classes }, el));
    }
    Split(panes, {
        direction: dir,
        gutterSize: 5,
        cursor: "pointer",
        sizes: config.sizes,
        onDragEnd: () => window.dispatchEvent(new Event("resize")),
    });
    return panes;
}

function createMenuIcon(icon, config, parent) {
    const attrs = {
        role: "tooltip", //needs to be tooltip for microtip to work.
    };
    const t = config.tooltip;
    if (t) {
        Object.assign(attrs, {
            "aria-label": t.text,
            "data-microtip-size": t.size || "small",
            "data-microtip-position": t.position || "bottom-left",
        });
    }

    const sp = createEl("span", attrs); //would changing to buttton make ir a11y compliant

    createEl(
        "i",
        {
            classes: ["ciview-menu-icon"].concat(icon.split(" ")),
            styles: {
                fontSize: config.size || "18px",
            },
        },
        sp,
    );
    if (config.func) {
        sp.addEventListener("click", (e) => config.func(e));
    }
    if (parent) {
        parent.append(sp);
    }
    return sp;
}
function addElProps(el, attrs) {
    for (const idx in attrs) {
        if ((idx === "styles" || idx === "style") && typeof attrs[idx] === "object") {
            for (const prop in attrs[idx]) {
                el.style[prop] = attrs[idx][prop];
            }
        } else if (idx === "text") {
            el.textContent = attrs[idx];
        } else if (idx === "classes") {
            for (const cl of attrs[idx]) {
                el.classList.add(cl);
            }
        } else {
            el.setAttribute(idx, attrs[idx]);
        }
    }
}

/**
 * Make a text input element that can be used to filter a list of items
 * @param {HTMLSelectElement} selectEl - the element to be filtered
 * @param {HTMLElement=} parent - (optional) the parent element to which the filter will be added
 * @returns {HTMLInputElement} - the filter element
 */
export function createFilterElement(selectEl, parent) {
    //consider using AutoComplete here
    const filter = createEl(
        "input",
        {
            placeholder: "Filter",
            type: "text",
            styles: {
                width: "4em",
                margin: "0.2em",
            },
        },
        parent,
    );
    filter.oninput = () => {
        const val = filter.value.toLowerCase().split(" ");
        for (const o of selectEl.options) {
            const filter = val.some((v) => o.text.toLowerCase().indexOf(v) === -1);
            if (filter) {
                o.style.display = "none";
            } else {
                o.style.display = "block";
            }
        }
    };
    return filter;
}

/**
 * Enables the children of an HTML element to be sorted by the user
 * via drag and drop
 * @param {HTMLElement} list - an HTML element whose immediate children will
 * be capable of being sorted.
 * @param {object} config
 * @param {string} [config.handle] - if each child element is to be only
 * dragged by a handle, then this should be a class name that the handle
 * contains.
 * @param {function} [config.sortEnded] - a function that is called when
 * a sorting action has finished. The function is supplied with the
 * current sorted list of child elements.
 */
function makeSortable(list, config = {}) {
    let dragged = null;
    let allowDrag = false;

    function handleDragOver(e) {
        //place the dragged item before the current item
        this.after(dragged);
        e.preventDefault();
        return false;
    }

    function handleDragEnd(e) {
        allowDrag = false;
        dragged.classList.remove("mdv-element-dragged");
        if (config.sortEnded) {
            config.sortEnded(Array.from(list.children).slice(1));
        }
    }

    function handleMouseDown(e) {
        //check whether item can be sorted
        if (!config.handle) {
            allowDrag = true;
        } else {
            if (e.target.classList.contains(config.handle)) {
                allowDrag = true;
            } else {
                allowDrag = false;
            }
        }
    }

    function handleDragStart(e) {
        if (!allowDrag) {
            e.preventDefault();
            return false;
        }
        dragged = this;
        this.classList.add("mdv-element-dragged");
        //blank drag image as the swapping of items
        //will give the appearance of dragging
        e.dataTransfer.setDragImage(createEl("div", {}), 0, 0);
    }

    //add dummy element, allows items to be dragged to the beginning
    //of the list i.e. after the dummy element
    const firstDummyEl = createEl("div", { styles: { height: "5px" } });
    list.prepend(firstDummyEl);
    firstDummyEl.addEventListener("dragover", handleDragOver);
    //add drag events to all child elements
    for (let n = 0; n < list.children.length; n++) {
        const item = list.children[n];
        item.draggable = true;
        item.addEventListener("dragover", handleDragOver);
        item.addEventListener("dragstart", handleDragStart);
        item.addEventListener("dragend", handleDragEnd);
        item.addEventListener("mousedown", handleMouseDown);
    }
}

function makeResizable(el, config = {}) {
    //already resizable
    if (el.__resizeinfo__) {
        return;
    }
    const ri = {
        resize: el.style.resize,
        overflow: el.style.overflow,
    };
    //document can change if in another window
    el.__doc__ = config.doc || document;
    // el.style.resize="both"; //standard resizer is sometimes visible when it shouldn't be.
    el.style.overflow = "visible";
    //el.style.zIndex="0";
    if (config.onresizeend) {
        ri.onresize = addResizeListener(
            el,
            (x, y) => {
                config.onresizeend(x, y);
            },
            config.onResizeStart,
        );
    }
    //workaround for Safari bug #50
    //https://codepen.io/jkasun/pen/QrLjXP
    //(actually, mostly copilot filling in very similar code...)

    // List of all the resizers based on directions and corresponding css classes
    const directions = [
        { dir: "n", className: "resizer-n" },
        { dir: "s", className: "resizer-s" },
        { dir: "e", className: "resizer-e" },
        { dir: "w", className: "resizer-w" },
        { dir: "ne", className: "resizer-ne" },
        { dir: "nw", className: "resizer-nw" },
        { dir: "se", className: "resizer-se" },
        { dir: "sw", className: "resizer-sw" },
    ];

    const resizeEls = directions.map(({ className }) => {
        const resizeEl = createEl(
            "div",
            {
                classes: [className],
            },
            el,
        );
        resizeEl.addEventListener("mousedown", initDrag, false);
        return resizeEl;
    });

    ri.removeChildren = () => {
        resizeEls.forEach((el) => {
            el.remove();
        });
    };
    function initDrag(e) {
        e.preventDefault();
        if (config.onResizeStart) config.onResizeStart();
        const startX = e.clientX;
        const startY = e.clientY;
        const startWidth = el.offsetWidth;
        const startHeight = el.offsetHeight;
        const startLeft = el.offsetLeft;
        const startTop = el.offsetTop;

        // Determine the direction from the target's class list
        const target = e.target;
        const directionObj = directions.find(({ className }) => target.classList.contains(className));
        const dir = directionObj ? directionObj.dir : "";

        function doDrag(e) {
            // Change in coordinates
            const dx = e.clientX - startX;
            const dy = e.clientY - startY;

            let newWidth = startWidth;
            let newHeight = startHeight;
            let newLeft = startLeft;
            let newTop = startTop;

            if (dir.includes("e")) {
                newWidth = startWidth + dx;
            }
            if (dir.includes("w")) {
                newWidth = startWidth - dx;
                newLeft = startLeft + dx;
            }
            if (dir.includes("s")) {
                newHeight = startHeight + dy;
            }
            if (dir.includes("n")) {
                newHeight = startHeight - dy;
                newTop = startTop + dy;
            }
            el.style.width = `${newWidth}px`;
            el.style.height = `${newHeight}px`;
            el.style.left = `${newLeft}px`;
            el.style.top = `${newTop}px`;
        }

        // Cleanup
        function stopDrag() {
            el.__doc__.documentElement.removeEventListener("mousemove", doDrag, false);
            el.__doc__.documentElement.removeEventListener("mouseup", stopDrag, false);
            if (config.onresizeend) {
                config.onresizeend(el.offsetWidth, el.offsetHeight, el.offsetLeft, el.offsetTop);
            }
        }

        el.__doc__.documentElement.addEventListener("mousemove", doDrag, false);
        el.__doc__.documentElement.addEventListener("mouseup", stopDrag, false);
    }

    el.__resizeinfo__ = ri;
}
function removeResizable(el) {
    if (!el.__resizeinfo__) {
        return;
    }
    const ri = el.__resizeinfo__;
    el.style.resize = ri.resize;
    el.style.overflow = ri.overflow;
    if (ri.onresize) {
        el.removeEventListener("mouseup", ri.onresize);
    }
    ri.removeChildren?.();
    el.__resizeinfo__ = undefined;
}

function removeDraggable(el) {
    const d = el.__draginfo__;
    if (!d) {
        return;
    }
    d.handle.style.cursor = d.cursor;
    el.style.position = d.position;
    d.handle.onmousedown = null;
    el.__draginfo__ = undefined;
}

class MDVProgress {
    constructor(el, config) {
        this.el = el;
        this.max = config.max || 100;
        el.classList.add("mdv-progress-outer");
        this.bar = createEl(
            "div",
            {
                classes: ["mdv-progress-inner"],
            },
            el,
        );
        this.text = createEl(
            "div",
            {
                classes: ["mdv-progress-text"],
            },
            el,
        );
        this.setText(config.text || "");
        this.setValue(config.value || 0);
    }

    setText(text) {
        this.text.textContent = text;
    }
    setValue(val) {
        this.bar.style.width = `${Math.round((val / this.max) * 100)}%`;
    }
}

/**
 * Make an element draggable
 * @param {HTMLElement} el - the element to be made draggable
 * @param {object} config
 * @param {HTMLElement} [config.handle] - the element that will be used to drag the main element
 * @param {HTMLElement} [config.contain] - the element that will contain the draggable element
 * @param {string} [config.y_axis] - if set, the element will only be draggable on the y-axis
 * @param {function} [config.ondragstart] - a function that is called when dragging starts
 * @param {function} [config.ondragend] - a function that is called when dragging ends
 * @param {Document} [config.doc] - the document object to which the event listeners
 * will be added. This is useful when the draggable element is in a different window
 * @param {boolean} [config.snapback] - if set, if the drag is released with the element off 
 * the side of its container, it will snap back to the edge. This has known issues and is
 * currently only being used for Dialogs.
 * @returns {void}
 */
function makeDraggable(el, config = {}) {
    if (!config.doc) {
        config.doc = document;
    }

    let pos1 = 0;
    let pos2 = 0;
    let pos3 = 0;
    let pos4 = 0;
    const handle = config.handle ? el.querySelector(config.handle) : el;
    let cont = null;
    if (config.contain) {
        cont = {
            dir: config.contain,
        };
    }

    handle.onmousedown = dragMouseDown;

    el.__draginfo__ = {
        handle: handle,
        cursor: handle.style.cursor,
        position: el.style.position,
    };
    handle.style.cursor = "move";
    el.style.position = "absolute";
    el.style.margin = "0px 0px 0px 0px";
    el.__doc__ = config.doc;

    function dragMouseDown(e) {
        e = e || window.event;
        e.preventDefault();

        //! current fix for the dropdown issue while dragging, check for any side effects
        if (document.activeElement) {
            document.activeElement.blur();
        }

        // get the mouse cursor position at startup:
        pos3 = e.clientX;
        pos4 = e.clientY;
        if (config.ondragstart) {
            config.ondragstart(el);
        }
        if (cont) {
            cont.p_bb = el.parentElement.getBoundingClientRect();
            cont.c_bb = el.getBoundingClientRect();
        }
        el.__doc__.onmouseup = closeDragElement;
        el.__doc__.onmousemove = elementDrag;
    }

    function elementDrag(e) {
        e = e || window.event;
        e.preventDefault();
        e.stopPropagation();
        // calculate the new cursor position:
        pos1 = pos3 - e.clientX;
        pos2 = pos4 - e.clientY;
        pos3 = e.clientX;
        pos4 = e.clientY;
        // set the element's new position:
        const nt = el.offsetTop - pos2;
        const nl = el.offsetLeft - pos1;
        if (cont) {
            if (nt < 0 || (nt + cont.c_bb.height > cont.p_bb.height && cont.dir !== "topleft")) {
                return;
            }
            if (nl < 0 || (nl + cont.c_bb.width > cont.p_bb.width && cont.dir !== "topleft")) {
                return;
            }
        }

        if (!config.y_axis) {
            el.style.top = `${nt}px`;
        }

        el.style.left = `${nl}px`;
    }

    function closeDragElement() {
        // stop moving when mouse button is released:
        if (config.ondragend) {
            config.ondragend();
        }

        if (config.snapback) {
            // if it is (nearly) off the side of the window, keep it within 100px
            const { top, right, left } = el.getBoundingClientRect();
            //! current logic for comparisons with parentBB is bad when applied to e.g. charts
            //need to think differently about relative coordinates
            const parentBB = el.parentElement.getBoundingClientRect();
            if (right < 100) {
                el.style.left = `${100 - el.offsetWidth}px`;
            }
            if (top > (parentBB.height - 100)) {
                el.style.top = `${parentBB.height - 100}px`;
            }
            if (left > (parentBB.width - 100)) {
                el.style.left = `${parentBB.width - 100}px`;
            }
            // if it is going out of the window at the top, bring the top to 0
            if (top < 0) {
                el.style.top = '0px';
            }
        }

        el.__doc__.onmouseup = null;
        el.__doc__.onmousemove = null;
    }
}

function addResizeListener(element, endCallback, startCallback) {
    const box = element.getBoundingClientRect();
    let width = box.width;
    let height = box.height;
    const list = (e) => {
        const box = element.getBoundingClientRect();
        if (box.width !== width || box.height !== height) {
            endCallback(box.width, box.height);
        }
        width = box.width;
        height = box.height;
    };
    if (startCallback) element.addEventListener("mousedown", startCallback);
    element.addEventListener("mouseup", list);
    return list;
}

function getElDim(el) {
    const rect = el.getBoundingClientRect();
    const scrollLeft = window.pageXOffset || document.documentElement.scrollLeft;
    const scrollTop = window.pageYOffset || document.documentElement.scrollTop;
    return {
        top: rect.top + scrollTop,
        left: rect.left + scrollLeft,
        height: rect.height,
        width: rect.width,
    };
}

export {
    createEl,
    createSVGEl,
    addResizeListener,
    makeDraggable,
    addElProps,
    makeResizable,
    removeDraggable,
    removeResizable,
    createMenuIcon,
    makeSortable,
    splitPane,
    getElDim,
    MDVProgress,
};
