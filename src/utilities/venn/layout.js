
import {intersectionArea, circleOverlap, circleCircleIntersection, distance} from './circleintersection.js';

function zeros(x){
	const r = new Array(x);
	for (let i = 0; i < x; ++i) {
		r[i] = 0; 
    }
    return r;
}
function zerosM(x,y){
	return zeros(x).map(() => zeros(y));
}

function norm2(a)  {
    return Math.sqrt(dot(a, a));
}

function scale(ret, value, c) {
    for (let i = 0; i < value.length; ++i) {
        ret[i] = value[i] * c;
    }
}
function weightedSum(ret, w1, v1, w2, v2) {
    for (let j = 0; j < ret.length; ++j) {
        ret[j] = w1 * v1[j] + w2 * v2[j];
    }
}

function gemv(output, A, x) {
    for (let i = 0; i < output.length; ++i) {
        output[i] = dot(A[i], x);
    }
}
function dot(a, b) {
    let ret = 0;
    for (let i = 0; i < a.length; ++i) {
        ret += a[i] * b[i];
    }
    return ret;
}


function wolfeLineSearch(f, pk, current, next, a, c1, c2) {
    // biome-ignore lint/style/useSingleVarDeclarator: don't mind too much
    let phi0 = current.fx, phiPrime0 = dot(current.fxprime, pk),
        phi = phi0, phi_old = phi0,
        phiPrime = phiPrime0,
        a0 = 0;

    a = a || 1;
    c1 = c1 || 1e-6;
    c2 = c2 || 0.1;

    function zoom(a_lo, a_high, phi_lo) {
        for (let iteration = 0; iteration < 16; ++iteration) {
            a = (a_lo + a_high)/2;
            weightedSum(next.x, 1.0, current.x, a, pk);
            phi = next.fx = f(next.x, next.fxprime);
            phiPrime = dot(next.fxprime, pk);

            if ((phi > (phi0 + c1 * a * phiPrime0)) ||
                (phi >= phi_lo)) {
                a_high = a;

            } else  {
                if (Math.abs(phiPrime) <= -c2 * phiPrime0) {
                    return a;
                }

                if (phiPrime * (a_high - a_lo) >=0) {
                    a_high = a_lo;
                }

                a_lo = a;
                phi_lo = phi;
            }
        }

        return 0;
    }

    for (let iteration = 0; iteration < 10; ++iteration) {
        weightedSum(next.x, 1.0, current.x, a, pk);
        phi = next.fx = f(next.x, next.fxprime);
        phiPrime = dot(next.fxprime, pk);
        if ((phi > (phi0 + c1 * a * phiPrime0)) ||
            (iteration && (phi >= phi_old))) {
            return zoom(a0, a, phi_old);
        }

        if (Math.abs(phiPrime) <= -c2 * phiPrime0) {
            return a;
        }

        if (phiPrime >= 0 ) {
            return zoom(a, a0, phi);
        }

        phi_old = phi;
        a0 = a;
        a *= 2;
    }

    return a;
}


function conjugateGradient(f, initial, params) {
    // allocate all memory up front here, keep out of the loop for perfomance
    // reasons
    // biome-ignore lint/style/useSingleVarDeclarator: this is a bit ugly but never mind
    let current = {x: initial.slice(), fx: 0, fxprime: initial.slice()},
        next = {x: initial.slice(), fx: 0, fxprime: initial.slice()},
        yk = initial.slice(),
        pk, temp,
        a = 1,
        maxIterations;

    params = params || {};
    maxIterations = params.maxIterations || initial.length * 20;

    current.fx = f(current.x, current.fxprime);
    pk = current.fxprime.slice();
    scale(pk, current.fxprime,-1);

    for (let i = 0; i < maxIterations; ++i) {
        a = wolfeLineSearch(f, pk, current, next, a);

        // todo: history in wrong spot?
        if (params.history) {
            params.history.push({x: current.x.slice(),
                                 fx: current.fx,
                                 fxprime: current.fxprime.slice(),
                                 alpha: a});
        }

        if (!a) {
            // faiiled to find point that satifies wolfe conditions.
            // reset direction for next iteration
            scale(pk, current.fxprime, -1);

        } else {
            // update direction using Polak�Ribiere CG method
            weightedSum(yk, 1, next.fxprime, -1, current.fxprime);

            const delta_k = dot(current.fxprime, current.fxprime);
            const beta_k = Math.max(0, dot(yk, next.fxprime) / delta_k);

            weightedSum(pk, beta_k, pk, -1, next.fxprime);

            temp = current;
            current = next;
            next = temp;
        }

        if (norm2(current.fxprime) <= 1e-5) {
            break;
        }
    }

    if (params.history) {
        params.history.push({x: current.x.slice(),
                             fx: current.fx,
                             fxprime: current.fxprime.slice(),
                             alpha: a});
    }

    return current;
}

function bisect(f, a, b, parameters) {
    parameters = parameters || {};
    // biome-ignore lint/style/useSingleVarDeclarator: <explanation>
    let maxIterations = parameters.maxIterations || 100,
        tolerance = parameters.tolerance || 1e-10,
        fA = f(a),
        fB = f(b),
        delta = b - a;

    if (fA * fB > 0) {
        throw "Initial bisect points must have opposite signs";
    }

    if (fA === 0) return a;
    if (fB === 0) return b;

    for (let i = 0; i < maxIterations; ++i) {
        delta /= 2;
        const mid = a + delta;
        const fMid = f(mid);

        if (fMid * fA >= 0) {
            a = mid;
        }

        if ((Math.abs(delta) < tolerance) || (fMid === 0)) {
            return mid;
        }
    }
    return a + delta;
}

function nelderMead(f, x0, parameters) {
    parameters = parameters || {};

    // biome-ignore lint/style/useSingleVarDeclarator: <explanation>
    let maxIterations = parameters.maxIterations || x0.length * 200,
        nonZeroDelta = parameters.nonZeroDelta || 1.05,
        zeroDelta = parameters.zeroDelta || 0.001,
        minErrorDelta = parameters.minErrorDelta || 1e-6,
        minTolerance = parameters.minErrorDelta || 1e-5,
        rho = (parameters.rho !== undefined) ? parameters.rho : 1,
        chi = (parameters.chi !== undefined) ? parameters.chi : 2,
        psi = (parameters.psi !== undefined) ? parameters.psi : -0.5,
        sigma = (parameters.sigma !== undefined) ? parameters.sigma : 0.5,
        maxDiff;

    // initialize simplex.
    const N = x0.length;
    const simplex = new Array(N + 1);
    simplex[0] = x0;
    simplex[0].fx = f(x0);
    simplex[0].id = 0;
    for (let i = 0; i < N; ++i) {
        const point = x0.slice();
        point[i] = point[i] ? point[i] * nonZeroDelta : zeroDelta;
        simplex[i+1] = point;
        simplex[i+1].fx = f(point);
        simplex[i+1].id = i+1;
    }

    function updateSimplex(value) {
        for (let i = 0; i < value.length; i++) {
            simplex[N][i] = value[i];
        }
        simplex[N].fx = value.fx;
    }

    const sortOrder = (a, b) => a.fx - b.fx;

    const centroid = x0.slice();
    const reflected = x0.slice();
    const contracted = x0.slice();
    const expanded = x0.slice();

    for (let iteration = 0; iteration < maxIterations; ++iteration) {
        simplex.sort(sortOrder);

        if (parameters.history) {
            // copy the simplex (since later iterations will mutate) and
            // sort it to have a consistent order between iterations
            const sortedSimplex = simplex.map((x) => {
                const state = x.slice();
                state.fx = x.fx;
                state.id = x.id;
                return state;
            });
            sortedSimplex.sort((a,b) => a.id - b.id);

            parameters.history.push({x: simplex[0].slice(),
                                     fx: simplex[0].fx,
                                     simplex: sortedSimplex});
        }

        maxDiff = 0;
        for (let i = 0; i < N; ++i) {
            maxDiff = Math.max(maxDiff, Math.abs(simplex[0][i] - simplex[1][i]));
        }

        if ((Math.abs(simplex[0].fx - simplex[N].fx) < minErrorDelta) &&
            (maxDiff < minTolerance)) {
            break;
        }

        // compute the centroid of all but the worst point in the simplex
        for (let i = 0; i < N; ++i) {
            centroid[i] = 0;
            for (let j = 0; j < N; ++j) {
                centroid[i] += simplex[j][i];
            }
            centroid[i] /= N;
        }

        // reflect the worst point past the centroid  and compute loss at reflected
        // point
        const worst = simplex[N];
        weightedSum(reflected, 1+rho, centroid, -rho, worst);
        reflected.fx = f(reflected);

        // if the reflected point is the best seen, then possibly expand
        if (reflected.fx < simplex[0].fx) {
            weightedSum(expanded, 1+chi, centroid, -chi, worst);
            expanded.fx = f(expanded);
            if (expanded.fx < reflected.fx) {
                updateSimplex(expanded);
            }  else {
                updateSimplex(reflected);
            }
        }

        // if the reflected point is worse than the second worst, we need to
        // contract
        else if (reflected.fx >= simplex[N-1].fx) {
            let shouldReduce = false;

            if (reflected.fx > worst.fx) {
                // do an inside contraction
                weightedSum(contracted, 1+psi, centroid, -psi, worst);
                contracted.fx = f(contracted);
                if (contracted.fx < worst.fx) {
                    updateSimplex(contracted);
                } else {
                    shouldReduce = true;
                }
            } else {
                // do an outside contraction
                weightedSum(contracted, 1-psi * rho, centroid, psi*rho, worst);
                contracted.fx = f(contracted);
                if (contracted.fx < reflected.fx) {
                    updateSimplex(contracted);
                } else {
                    shouldReduce = true;
                }
            }

            if (shouldReduce) {
                // if we don't contract here, we're done
                if (sigma >= 1) break;

                // do a reduction
                for (i = 1; i < simplex.length; ++i) {
                    weightedSum(simplex[i], 1 - sigma, simplex[0], sigma, simplex[i]);
                    simplex[i].fx = f(simplex[i]);
                }
            }
        } else {
            updateSimplex(reflected);
        }
    }

    simplex.sort(sortOrder);
    return {fx : simplex[0].fx,
            x : simplex[0]};
}


/** given a list of set objects, and their corresponding overlaps.
updates the (x, y, radius) attribute on each set such that their positions
roughly correspond to the desired overlaps */
function venn(areas, parameters) {
    parameters = parameters || {};
    parameters.maxIterations = parameters.maxIterations || 500;
    const initialLayout = parameters.initialLayout || bestInitialLayout;
    const loss = parameters.lossFunction || lossFunction;

    // add in missing pairwise areas as having 0 size
    areas = addMissingAreas(areas);

    // initial layout is done greedily
    const circles = initialLayout(areas, parameters);

    // transform x/y coordinates to a vector to optimize
    const initial = [];
    const setids = [];
    let setid;
    for (setid in circles) {
        if (circles.hasOwnProperty(setid)) {
            initial.push(circles[setid].x);
            initial.push(circles[setid].y);
            setids.push(setid);
        }
    }

    // optimize initial layout from our loss function
    let totalFunctionCalls = 0;
    const solution = nelderMead(
        (values) => {
            totalFunctionCalls += 1;
            const current = {};
            for (let i = 0; i < setids.length; ++i) {
                const setid = setids[i];
                current[setid] = {x: values[2 * i],
                                  y: values[2 * i + 1],
                                  radius : circles[setid].radius,
                                 // size : circles[setid].size
                                 };
            }
            return loss(current, areas);
        },
        initial,
        parameters);

    // transform solution vector back to x/y points
    const positions = solution.x;
    for (let i = 0; i < setids.length; ++i) {
        setid = setids[i];
        circles[setid].x = positions[2 * i];
        circles[setid].y = positions[2 * i + 1];
    }

    return circles;
}

const SMALL = 1e-10;

/** Returns the distance necessary for two circles of radius r1 + r2 to
have the overlap area 'overlap' */
function distanceFromIntersectArea(r1, r2, overlap) {
    // handle complete overlapped circles
    if (Math.min(r1, r2) * Math.min(r1,r2) * Math.PI <= overlap + SMALL) {
        return Math.abs(r1 - r2);
    }

    return bisect((distance) => circleOverlap(r1, r2, distance) - overlap, 0, r1 + r2);
}

/** Missing pair-wise intersection area data can cause problems:
 treating as an unknown means that sets will be laid out overlapping,
 which isn't what people expect. To reflect that we want disjoint sets
 here, set the overlap to 0 for all missing pairwise set intersections */
function addMissingAreas(areas) {
    areas = areas.slice();

    // two circle intersections that aren't defined
    const ids = [];
    const pairs = {};
    let i;
    let j;
    let a;
    let b;
    for (i = 0; i < areas.length; ++i) {
        const area = areas[i];
        if (area.sets.length === 1) {
            ids.push(area.sets[0]);
        } else if (area.sets.length === 2) {
            a = area.sets[0];
            b = area.sets[1];
            pairs[[a, b]] = true;
            pairs[[b, a]] = true;
        }
    }
    ids.sort((a, b) => a > b);

    for (i = 0; i < ids.length; ++i) {
        a = ids[i];
        for (j = i + 1; j < ids.length; ++j) {
            b = ids[j];
            if (!([a, b] in pairs)) {
                areas.push({'sets': [a, b],
                            'size': 0});
            }
        }
    }
    return areas;
}

/// Returns two matrices, one of the euclidean distances between the sets
/// and the other indicating if there are subset or disjoint set relationships
function getDistanceMatrices(areas, sets, setids) {
    // initialize an empty distance matrix between all the points
    const distances = zerosM(sets.length, sets.length);
    const constraints = zerosM(sets.length, sets.length);

    // compute required distances between all the sets such that
    // the areas match
    areas.filter((x) => x.sets.length === 2)
        .map((current) => {
        const left = setids[current.sets[0]];
        const right = setids[current.sets[1]];
        const r1 = Math.sqrt(sets[left].size / Math.PI);
        const r2 = Math.sqrt(sets[right].size / Math.PI);
        const distance = distanceFromIntersectArea(r1, r2, current.size);

        distances[left][right] = distances[right][left] = distance;

        // also update constraints to indicate if its a subset or disjoint
        // relationship
        let c = 0;
        if (current.size + 1e-10 >= Math.min(sets[left].size,
                                             sets[right].size)) {
            c = 1;
        } else if (current.size <= 1e-10) {
            c = -1;
        }
        constraints[left][right] = constraints[right][left] = c;
    });

    return {distances: distances, constraints: constraints};
}

/// computes the gradient and loss simulatenously for our constrained MDS optimizer
function constrainedMDSGradient(x, fxprime, distances, constraints) {
    let loss = 0;
    let i;
    for (i = 0; i < fxprime.length; ++i) {
        fxprime[i] = 0;
    }

    for (i = 0; i < distances.length; ++i) {
        const xi = x[2 * i];
        const yi = x[2 * i + 1];
        for (let j = i + 1; j < distances.length; ++j) {
            const xj = x[2 * j];
            const yj = x[2 * j + 1];
            const dij = distances[i][j];
            const constraint = constraints[i][j];

            const squaredDistance = (xj - xi) * (xj - xi) + (yj - yi) * (yj - yi);
            const distance = Math.sqrt(squaredDistance);
            const delta = squaredDistance - dij * dij;

            if (((constraint > 0) && (distance <= dij)) ||
                ((constraint < 0) && (distance >= dij))) {
                continue;
            }

            loss += 2 * delta * delta;

            fxprime[2*i]     += 4 * delta * (xi - xj);
            fxprime[2*i + 1] += 4 * delta * (yi - yj);

            fxprime[2*j]     += 4 * delta * (xj - xi);
            fxprime[2*j + 1] += 4 * delta * (yj - yi);
        }
    }
    return loss;
}

/// takes the best working variant of either constrained MDS or greedy
export function bestInitialLayout(areas, params) {
    let initial = greedyLayout(areas, params);
    const loss = params.lossFunction || lossFunction;

    // greedylayout is sufficient for all 2/3 circle cases. try out
    // constrained MDS for higher order problems, take its output
    // if it outperforms. (greedy is aesthetically better on 2/3 circles
    // since it axis aligns)
    if (areas.length >= 8) {
        const constrained  = constrainedMDSLayout(areas, params);
        const constrainedLoss = loss(constrained, areas);
        const greedyLoss = loss(initial, areas);

        if (constrainedLoss + 1e-8 < greedyLoss) {
            initial = constrained;
        }
    }
    return initial;
}

/// use the constrained MDS variant to generate an initial layout
function constrainedMDSLayout(areas, params) {
    params = params || {};
    const restarts = params.restarts || 10;

    // bidirectionally map sets to a rowid  (so we can create a matrix)
    const sets = [];
    const setids = {};
    let i;
    for (i = 0; i < areas.length; ++i ) {
        const area = areas[i];
        if (area.sets.length === 1) {
            setids[area.sets[0]] = sets.length;
            sets.push(area);
        }
    }

    const matrices = getDistanceMatrices(areas, sets, setids);
    let distances = matrices.distances;
    const constraints = matrices.constraints;

    // keep distances bounded, things get messed up otherwise.
    // TODO: proper preconditioner?
    const norm = norm2(distances.map(norm2))/(distances.length);
    distances = distances.map((row) => row.map((value) => value / norm));

    const obj = (x, fxprime) => constrainedMDSGradient(x, fxprime, distances, constraints);

    let best;
    let current;
    for (i = 0; i < restarts; ++i) {
        const initial = zeros(distances.length*2).map(Math.random);

        current = conjugateGradient(obj, initial, params);
        if (!best || (current.fx < best.fx)) {
            best = current;
        }
    }
    const positions = best.x;

    // translate rows back to (x,y,radius) coordinates
    const circles = {};
    for (i = 0; i < sets.length; ++i) {
        const set = sets[i];
        circles[set.sets[0]] = {
            x: positions[2*i] * norm,
            y: positions[2*i + 1] * norm,
            radius:  Math.sqrt(set.size / Math.PI)
        };
    }

    if (params.history) {
        for (i = 0; i < params.history.length; ++i) {
            scale(params.history[i].x, norm);
        }
    }
    return circles;
}

/** Lays out a Venn diagram greedily, going from most overlapped sets to
least overlapped, attempting to position each new set such that the
overlapping areas to already positioned sets are basically right */
function greedyLayout(areas, params) {
    const loss = params?.lossFunction ? params.lossFunction : lossFunction;
    // define a circle for each set
    const circles = {};
    const setOverlaps = {};
    let set;
    for (let i = 0; i < areas.length; ++i) {
        const area = areas[i];
        if (area.sets.length === 1) {
            set = area.sets[0];
            circles[set] = {x: 1e10, y: 1e10,
                            rowid: circles.length,
                            size: area.size,
                            radius: Math.sqrt(area.size / Math.PI)};
            setOverlaps[set] = [];
        }
    }
    areas = areas.filter((a) => a.sets.length === 2);

    // map each set to a list of all the other sets that overlap it
    for (i = 0; i < areas.length; ++i) {
        const current = areas[i];
        let weight = current.hasOwnProperty('weight') ? current.weight : 1.0;
        const left = current.sets[0];
        const right = current.sets[1];

        // completely overlapped circles shouldn't be positioned early here
        if (current.size + SMALL >= Math.min(circles[left].size,
                                             circles[right].size)) {
            weight = 0;
        }

        setOverlaps[left].push ({set:right, size:current.size, weight:weight});
        setOverlaps[right].push({set:left,  size:current.size, weight:weight});
    }

    // get list of most overlapped sets
    const mostOverlapped = [];
    for (set in setOverlaps) {
        if (setOverlaps.hasOwnProperty(set)) {
            let size = 0;
            for (i = 0; i < setOverlaps[set].length; ++i) {
                size += setOverlaps[set][i].size * setOverlaps[set][i].weight;
            }

            mostOverlapped.push({set: set, size:size});
        }
    }

    // sort by size desc
    function sortOrder(a,b) {
        return b.size - a.size;
    }
    mostOverlapped.sort(sortOrder);

    // keep track of what sets have been laid out
    const positioned = {};
    function isPositioned(element) {
        return element.set in positioned;
    }

    // adds a point to the output
    function positionSet(point, index) {
        circles[index].x = point.x;
        circles[index].y = point.y;
        positioned[index] = true;
    }

    // add most overlapped set at (0,0)
    positionSet({x: 0, y: 0}, mostOverlapped[0].set);

    // get distances between all points. TODO, necessary?
    // answer: probably not
    // var distances = venn.getDistanceMatrices(circles, areas).distances;
    for (i = 1; i < mostOverlapped.length; ++i) {
        const setIndex = mostOverlapped[i].set;
        const overlap = setOverlaps[setIndex].filter(isPositioned);
        set = circles[setIndex];
        overlap.sort(sortOrder);

        if (overlap.length === 0) {
            // this shouldn't happen anymore with addMissingAreas
            throw "ERROR: missing pairwise overlap information";
        }

        const points = [];
        for (let j = 0; j < overlap.length; ++j) {
            // get appropriate distance from most overlapped already added set
            const p1 = circles[overlap[j].set];
            const d1 = distanceFromIntersectArea(set.radius, p1.radius,
                                               overlap[j].size);

            // sample positions at 90 degrees for maximum aesthetics
            points.push({x : p1.x + d1, y : p1.y});
            points.push({x : p1.x - d1, y : p1.y});
            points.push({y : p1.y + d1, x : p1.x});
            points.push({y : p1.y - d1, x : p1.x});

            // if we have at least 2 overlaps, then figure out where the
            // set should be positioned analytically and try those too
            for (let k = j + 1; k < overlap.length; ++k) {
                const p2 = circles[overlap[k].set];
                const d2 = distanceFromIntersectArea(set.radius, p2.radius,
                                                   overlap[k].size);

                const extraPoints = circleCircleIntersection(
                    { x: p1.x, y: p1.y, radius: d1},
                    { x: p2.x, y: p2.y, radius: d2});

                for (let l = 0; l < extraPoints.length; ++l) {
                    points.push(extraPoints[l]);
                }
            }
        }

        // we have some candidate positions for the set, examine loss
        // at each position to figure out where to put it at
        let bestLoss = 1e50;
        let bestPoint = points[0];
        for (j = 0; j < points.length; ++j) {
            circles[setIndex].x = points[j].x;
            circles[setIndex].y = points[j].y;
            const localLoss = loss(circles, areas);
            if (localLoss < bestLoss) {
                bestLoss = localLoss;
                bestPoint = points[j];
            }
        }

        positionSet(bestPoint, setIndex);
    }

    return circles;
}

/** Given a bunch of sets, and the desired overlaps between these sets - computes
the distance from the actual overlaps to the desired overlaps. Note that
this method ignores overlaps of more than 2 circles */
function lossFunction(sets, overlaps) {
    let output = 0;

    function getCircles(indices) {
        return indices.map((i) => sets[i]);
    }

    for (let i = 0; i < overlaps.length; ++i) {
        const area = overlaps[i];
        let overlap;
        if (area.sets.length === 1) {
            continue;
        }if (area.sets.length === 2) {
            const left = sets[area.sets[0]];
            const right = sets[area.sets[1]];
            overlap = circleOverlap(left.radius, right.radius,
                                    distance(left, right));
        } else {
            overlap = intersectionArea(getCircles(area.sets));
        }

        const weight = area.hasOwnProperty('weight') ? area.weight : 1.0;
        output += weight * (overlap - area.size) * (overlap - area.size);
    }

    return output;
}

// orientates a bunch of circles to point in orientation
function orientateCircles(circles, orientation, orientationOrder) {
    if (orientationOrder === null) {
        circles.sort((a, b) => b.radius - a.radius);
    } else {
        circles.sort(orientationOrder);
    }

    let i;
    // shift circles so largest circle is at (0, 0)
    if (circles.length > 0) {
        const largestX = circles[0].x;
        const largestY = circles[0].y;

        for (i = 0; i < circles.length; ++i) {
            circles[i].x -= largestX;
            circles[i].y -= largestY;
        }
    }

    if (circles.length === 2) {
        // if the second circle is a subset of the first, arrange so that
        // it is off to one side. hack for https://github.com/benfred/venn.js/issues/120
        const dist = distance(circles[0], circles[1]);
        if (dist < Math.abs(circles[1].radius - circles[0].radius)) {
            circles[1].x = circles[0].x + circles[0].radius - circles[1].radius - 1e-10;
            circles[1].y = circles[0].y;
        }
    }

    // rotate circles so that second largest is at an angle of 'orientation'
    // from largest
    if (circles.length > 1) {
        const rotation = Math.atan2(circles[1].x, circles[1].y) - orientation;
        const c = Math.cos(rotation);
        const s = Math.sin(rotation);
        let x;
        let y;

        for (i = 0; i < circles.length; ++i) {
            x = circles[i].x;
            y = circles[i].y;
            circles[i].x = c * x - s * y;
            circles[i].y = s * x + c * y;
        }
    }

    // mirror solution if third solution is above plane specified by
    // first two circles
    if (circles.length > 2) {
        let angle = Math.atan2(circles[2].x, circles[2].y) - orientation;
        while (angle < 0) { angle += 2* Math.PI; }
        while (angle > 2*Math.PI) { angle -= 2* Math.PI; }
        if (angle > Math.PI) {
            const slope = circles[1].y / (1e-10 + circles[1].x);
            for (i = 0; i < circles.length; ++i) {
                const d = (circles[i].x + slope * circles[i].y) / (1 + slope*slope);
                circles[i].x = 2 * d - circles[i].x;
                circles[i].y = 2 * d * slope - circles[i].y;
            }
        }
    }
}

function disjointCluster(circles) {
    // union-find clustering to get disjoint sets
    circles.map((circle) => { circle.parent = circle; });

    // path compression step in union find
    function find(circle) {
        if (circle.parent !== circle) {
            circle.parent = find(circle.parent);
        }
        return circle.parent;
    }

    function union(x, y) {
        const xRoot = find(x);
        const yRoot = find(y);
        xRoot.parent = yRoot;
    }

    // get the union of all overlapping sets
    for (let i = 0; i < circles.length; ++i) {
        for (let j = i + 1; j < circles.length; ++j) {
            const maxDistance = circles[i].radius + circles[j].radius;
            if (distance(circles[i], circles[j]) + 1e-10 < maxDistance) {
                union(circles[j], circles[i]);
            }
        }
    }

    // find all the disjoint clusters and group them together
    const disjointClusters = {};
    let setid;
    for (i = 0; i < circles.length; ++i) {
        setid = find(circles[i]).parent.setid;
        if (!(setid in disjointClusters)) {
            disjointClusters[setid] = [];
        }
        disjointClusters[setid].push(circles[i]);
    }

    // cleanup bookkeeping
    circles.map((circle) => { circle.parent = undefined; });

    // return in more usable form
    const ret = [];
    for (setid in disjointClusters) {
        if (disjointClusters.hasOwnProperty(setid)) {
            ret.push(disjointClusters[setid]);
        }
    }
    return ret;
}

function getBoundingBox(circles) {
    const minMax = (d) => {
        const hi = Math.max.apply(null, circles.map(
                                (c) => c[d] + c.radius ));
        const lo = Math.min.apply(null, circles.map(
                                (c) => c[d] - c.radius ));
        return {max:hi, min:lo};
    };

    return {xRange: minMax('x'), yRange: minMax('y')};
}

function normalizeSolution(solution, orientation, orientationOrder) {
    if (orientation === null){
        orientation = Math.PI/2;
    }

    // work with a list instead of a dictionary, and take a copy so we
    // don't mutate input
    let circles = [];
    let i;
    let setid;
    for (setid in solution) {
        if (solution.hasOwnProperty(setid)) {
            const previous = solution[setid];
            circles.push({x: previous.x,
                          y: previous.y,
                          radius: previous.radius,
                          setid: setid});
        }
    }

    // get all the disjoint clusters
    const clusters = disjointCluster(circles);

    // orientate all disjoint sets, get sizes
    for (i = 0; i < clusters.length; ++i) {
        orientateCircles(clusters[i], orientation, orientationOrder);
        const bounds = getBoundingBox(clusters[i]);
        clusters[i].size = (bounds.xRange.max - bounds.xRange.min) * (bounds.yRange.max - bounds.yRange.min);
        clusters[i].bounds = bounds;
    }
    clusters.sort((a, b) => b.size - a.size);

    // orientate the largest at 0,0, and get the bounds
    circles = clusters[0];
    let returnBounds = circles.bounds;

    const spacing = (returnBounds.xRange.max - returnBounds.xRange.min)/50;

    function addCluster(cluster, right, bottom) {
        if (!cluster) return;

        const bounds = cluster.bounds;
        let xOffset;
        let yOffset;
        let centreing;

        if (right) {
            xOffset = returnBounds.xRange.max  - bounds.xRange.min + spacing;
        } else {
            xOffset = returnBounds.xRange.max  - bounds.xRange.max;
            centreing = (bounds.xRange.max - bounds.xRange.min) / 2 -
                        (returnBounds.xRange.max - returnBounds.xRange.min) / 2;
            if (centreing < 0) xOffset += centreing;
        }

        if (bottom) {
            yOffset = returnBounds.yRange.max  - bounds.yRange.min + spacing;
        } else {
            yOffset = returnBounds.yRange.max  - bounds.yRange.max;
            centreing = (bounds.yRange.max - bounds.yRange.min) / 2 -
                        (returnBounds.yRange.max - returnBounds.yRange.min) / 2;
            if (centreing < 0) yOffset += centreing;
        }

        for (let j = 0; j < cluster.length; ++j) {
            cluster[j].x += xOffset;
            cluster[j].y += yOffset;
            circles.push(cluster[j]);
        }
    }

    let index = 1;
    while (index < clusters.length) {
        addCluster(clusters[index], true, false);
        addCluster(clusters[index+1], false, true);
        addCluster(clusters[index+2], true, true);
        index += 3;

        // have one cluster (in top left). lay out next three relative
        // to it in a grid
        returnBounds = getBoundingBox(circles);
    }

    // convert back to solution form
    const ret = {};
    for (i = 0; i < circles.length; ++i) {
        ret[circles[i].setid] = circles[i];
    }
    return ret;
}

/** Scales a solution from venn.venn or venn.greedyLayout such that it fits in
a rectangle of width/height - with padding around the borders. also
centers the diagram in the available space at the same time */
function scaleSolution(solution, width, height, padding) {
    const circles = [];
    const setids = [];
    for (const setid in solution) {
        if (solution.hasOwnProperty(setid)) {
            setids.push(setid);
            circles.push(solution[setid]);
        }
    }

    width -= 2*padding;
    height -= 2*padding;

    const bounds = getBoundingBox(circles);
    const xRange = bounds.xRange;
    const yRange = bounds.yRange;

    if ((xRange.max === xRange.min) ||
        (yRange.max === yRange.min)) {
        console.log("not scaling solution: zero size detected");
        return solution;
    }

    const xScaling = width  / (xRange.max - xRange.min);
    const yScaling = height / (yRange.max - yRange.min);
    const scaling = Math.min(yScaling, xScaling);
    // while we're at it, center the diagram too
    const xOffset = (width -  (xRange.max - xRange.min) * scaling) / 2;
    const yOffset = (height - (yRange.max - yRange.min) * scaling) / 2;

    const scaled = {};
    for (let i = 0; i < circles.length; ++i) {
        const circle = circles[i];
        scaled[setids[i]] = {
            radius: scaling * circle.radius,
            x: padding + xOffset + (circle.x - xRange.min) * scaling,
            y: padding + yOffset + (circle.y - yRange.min) * scaling,
        };
    }

    return scaled;
}


export {venn, lossFunction, normalizeSolution, scaleSolution,nelderMead};