import { makeDraggable } from "./Elements";
class ImagePanZoom {
    constructor(container, image, config = {}) {
        this.container = container;

        this.img = new Image();

        this.img.style.position = "absolute";
        this.img.style.maxWidth= 'none'

        this.img.onload = (e) => {
            this.container.append(this.img);
            this.fit();
        };
        makeDraggable(this.img);
        this.container.style.overflow = "hidden";
        this.container.addEventListener("wheel", (event) => {
            event.preventDefault();
            this.zoom(Math.sign(event.deltaY) > 0 ? -1 : 1, event);
        });
        this.setImage(image);
    }

    setImage(url) {
        this.img.src = url;
    }

    fit() {
        const box = this.container.getBoundingClientRect();
        const ratio =  this.img.naturalWidth / this.img.naturalHeight;
        this.img.style.height = `${box.height}px`;
        this.img.style.width = `${this.img.height * ratio}px`;
        if (this.img.width > box.width) {
            this.img.style.width = `${box.width}px`;
            this.img.style.height = `${this.img.width / ratio}px`;
            this.img.style.top = `${(box.height - this.img.height) / 2}px`;
            this.img.style.left = "0px";
        } else {
            this.img.style.top = "0px";
            this.img.style.left = `${(box.width - this.img.width) / 2}px`;
        }
    }

    zoom(amount, event) {
        amount = amount > 0 ? 1.1 : 0.9;
        const cbox = this.container.getBoundingClientRect();
        const ibox = this.img.getBoundingClientRect();
        const dx =
            (event.clientX - cbox.left - this.img.offsetLeft) * (amount - 1);
        const dy =
            (event.clientY - cbox.top - this.img.offsetTop) * (amount - 1);

        this.img.style.width = `${ibox.width * amount}px`;
        this.img.style.height = `${ibox.height * amount}px`;
        this.img.style.left = `${this.img.offsetLeft - dx}px`;
        this.img.style.top = `${this.img.offsetTop - dy}px`;
    }
}

export default ImagePanZoom;
