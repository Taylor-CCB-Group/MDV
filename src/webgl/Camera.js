import * as glMatrix from "gl-matrix";

function damp(x) {
    const xd = x * 0.9;
    if (Math.abs(xd) < 0.1) {
        return 0;
    }
    return xd;
}

function clamp(x, lo, hi) {
    return Math.min(Math.max(x, lo), hi);
}

class Camera {
    constructor(props = {}) {
        this.cameraState = {
            view: glMatrix.mat4.identity(new Float32Array(16)),
            projection: glMatrix.mat4.identity(new Float32Array(16)),
            center: props.center || [0, 0, 0], //new Float32Array(props.center || 3),
            theta: props.theta || 0,
            phi: props.phi || 0,
            distance: Math.log(props.distance || 1000),
            eye: new Float32Array(3),
            up: new Float32Array(props.up || [0, 1, 0]),
        };
        this.right = new Float32Array([1, 0, 0]);
        this.front = new Float32Array([0, 0, 1]);

        this.minDistance = Math.log(
            "minDistance" in props ? props.minDistance : 0.5,
        );
        this.maxDistance = Math.log(
            "maxDistance" in props ? props.maxDistance : 100000000000000,
        );

        this.dtheta = 0;
        this.dphi = 0;
        this.ddistance = 0;

        this.updateCamera();
    }

    mouseChange(dx, dy) {
        const w = Math.max(this.cameraState.distance, 0.5);

        this.dtheta += w * dx;
        this.dphi += w * dy;
        this.updateCamera();
    }

    getProjection() {
        return {
            projection: this.cameraState.projection,
            view: this.cameraState.view,
            distance: this.cameraState.distance,
            eye: this.cameraState.eye,
        };
    }

    mouseWheel(dx, dy) {
        this.ddistance += dy / window.innerHeight;
        this.updateCamera();
    }

    updateCamera() {
        const center = this.cameraState.center;
        const eye = this.cameraState.eye;
        const up = this.cameraState.up;

        this.cameraState.theta += this.dtheta;
        this.cameraState.phi = clamp(
            this.cameraState.phi + this.dphi,
            -Math.PI / 2.0,
            Math.PI / 2.0,
        );
        this.cameraState.distance = clamp(
            this.cameraState.distance + this.ddistance,
            this.minDistance,
            this.maxDistance,
        );

        this.dtheta = damp(this.dtheta);
        this.dphi = damp(this.dphi);
        this.ddistance = damp(this.ddistance);

        const theta = this.cameraState.theta;
        const phi = this.cameraState.phi;
        //console.log(phi+":"+theta)
        const r = Math.exp(this.cameraState.distance);

        const vf = r * Math.sin(theta) * Math.cos(phi);
        const vr = r * Math.cos(theta) * Math.cos(phi);
        const vu = r * Math.sin(phi);

        for (let i = 0; i < 3; ++i) {
            eye[i] =
                center[i] +
                vf * this.front[i] +
                vr * this.right[i] +
                vu * up[i];
        }

        const test = [60 * Math.cos(0.01), 30, 60 * Math.sin(0.01)];

        glMatrix.mat4.lookAt(this.cameraState.view, eye, center, up);
    }
}

export { Camera };
