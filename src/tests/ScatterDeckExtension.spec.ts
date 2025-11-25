// hope to get back to this once deck test-utils config issues are figured out.

// import { OrthographicView } from '@deck.gl/core';
// import { expect, test, vi } from 'vitest';
// import { generateLayerTests, testLayer } from '@deck.gl/test-utils';
// import { ScatterSquareExtension, ScatterDensityExension } from '../webgl/ScatterDeckExtension';
// import { type ScatterplotLayerProps, ScatterplotLayer } from "@deck.gl/layers";

// interface Position {
//   x: number;
//   y: number;
//   z: number;
// }

// test('ScatterSquareExtension', (t) => {
//   const view = new OrthographicView({
//     id: 'test',
//     width: 4,
//     height: 4,
//     // target: [2, 2, 0],
//     // zoom: 0,
//     viewState: {
//       target: [2, 2, 0],
//       zoom: 0,
//     }
//   });
//   const testCases = generateLayerTests({
//     Layer: ScatterplotLayer<Position>,
//     assert: (value, _msg) => expect(value).toBeTruthy(),
//     sampleProps: {
//       extensions: [new ScatterSquareExtension()],
//       data: [
//         { x: 0, y: 0, z: 0 },
//         { x: 1, y: 1, z: 1 },
//       ],
//       getPosition: (d: { x: number, y: number, z: number }) => [d.x, d.y, d.z] as [number, number, number],
//     },
//     onAfterUpdate: () => {
//       t.annotate('scatterplot with ScatterSquareExtension onAfterUpdate');
//     },
//   });
//   const viewport = view.makeViewport({
//     width: 4,
//     height: 4,
//     viewState: {
//       target: [2, 2, 0],
//       zoom: 0,
//     }
//   });
//   expect(viewport).not.toBeNull();
//   if (!viewport) throw new Error('viewport is null');
//   testLayer({
//     Layer: ScatterplotLayer<Position>,
//     testCases,
//     onError: err => {
//       // if (err) console.error(err);
//       // expect(err).toBeFalsy();
//     },
//     viewport,
//   })
// })