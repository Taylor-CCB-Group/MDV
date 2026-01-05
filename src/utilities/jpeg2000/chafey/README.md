# OpenJpegJS for raw codestreams

The code in this folder is forked from https://github.com/chafey/openjpegjs under MIT license

> JS/WebAssembly build of OpenJPEG
>
> NOTE - a forked version of OpenJPEG is currently used which has some changes to allow partial bitstream decoding

https://github.com/xinaesthete/openjpegjs/tree/arm64_build

Forked again to update build, emit ES6 & types, and patch endianness.

Once we have things working, we may try to publish a version to npm in a more civilized manner.

## Current status:

Assumes that images are encoded as BigEndian - this is what will be output by `raw2ometiff`, which for the time-being is our canonical way of generating these images.

Endianness should be determined from the first two bytes of the file, but I can't quite see how to get that information out of `geotiff.js` as of now.