/**
 * Creates a context menu the - the function should return a list describing the menu
 * The list should contain objects with text,icon,ghosted,func and subitems
 */
class ContextMenu {
    /**
     * @param {function} func - The method that will be used to construct the menu,
     * based on the data passed to the show method
     */
    constructor(func) {
        this.setItemFunction = func;
        this.menus = [];
    }

    /**
     * Will show the actual menu
     * @param {Event} e - The event that triggers the context menu
     * required to position the menu and also the menu is attached to event target
     * @param {any} data - optional, the menu will be constructed based on the function
     * given to the construtor and the data passed in this method
     *
     */
    show(e, data) {
        if (this.menus.length > 0) {
            this.remove();
        }
        const items = this.setItemFunction(data);

        this._addMenu(items, e.pageX, e.pageY);
        console.log(e.clientY);
        //stop propagation before listener is added
        e.stopPropagation();
        //remove the menu if user clicks outside
        this.removeFunc = () => {
            this.remove();
        };
        this._getDocument().body.addEventListener("click", this.removeFunc);
    }

    /**
     * closes the menu and removes it from the dom will be called uf the user
     * clicks outside the menu
     */
    remove() {
        for (const m of this.menus) {
            m.remove();
        }
        this.menus = [];
        if (this.removeFunc) {
            this._getDocument().body.removeEventListener(
                "click",
                this.removeFunc,
            );
        }
    }

    _addMenu(items, left, top) {
        const menu = document.createElement("div");
        const doc = this._getDocument();

        const index = this.menus.length;
        menu.classList.add("ciview-ctm-main");
        for (const item of items) {
            this._addItem(menu, item, index);
        }
        doc.body.append(menu);
        const rect = menu.getBoundingClientRect();
        if (left + rect.width > doc.body.clientWidth) {
            left = doc.body.clientWidth - rect.width - 3;
        }
        menu.style.left = `${left}px`;
        menu.style.top = `${top}px`;
        this.menus.push(menu);
    }
    _getDocument() {
        return this.__doc__ ? this.__doc__ : document;
    }

    _showSubmenu(el, subitems, index) {
        //remove all open menus below this level
        for (let n = index + 1; n < this.menus.length; n++) {
            this.menus[n].remove();
        }
        this.menus = this.menus.slice(0, index + 1);

        //work out where to place submenu
        const box = el.getBoundingClientRect();
        this._addMenu(subitems, box.left + box.width + 1, box.top);
    }

    _addItem(menu, item, index) {
        const el = document.createElement("div");
        el.classList.add("ciview-ctm-item");
        if (item.ghosted) {
            el.classList.add("ciview-ctm-ghosted");
        } else {
            el.classList.add("ciview-ctm-active");
        }

        if (item.icon) {
            const icon = document.createElement("i");
            for (const i of item.icon.split(" ")) {
                icon.classList.add(i);
            }
            icon.classList.add("ciview-ctm-left-icon");
            el.appendChild(icon);
        }

        const sp = document.createElement("span");
        sp.appendChild(document.createTextNode(item.text));
        el.append(sp);
        if (item.subitems) {
            const caret = document.createElement("i");
            caret.classList.add(
                "fas",
                "fa-caret-right",
                "ciview-ctm-right-icon",
            );
            el.appendChild(caret);
        }
        el.addEventListener("click", (e) => {
            if (item.func && !item.ghosted) {
                item.func(e);
            }
            if (item.subitems) {
                this._showSubmenu(el, item.subitems, index);
                for (const i of menu.children) {
                    i.classList.remove("ciview-ctm-selected");
                }
                el.classList.add("ciview-ctm-selected");
            } else {
                this.remove();
            }
            e.stopPropagation();
        });
        menu.append(el);
    }
}

export { ContextMenu };
