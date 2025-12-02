//! Commented due to incompatibility with React 19
// /// <reference types="@welldone-software/why-did-you-render" />
// import React from "react";

// async function loadWdyr() {
//     console.log("Development mode...");
//     const { default: whyDidYouRender } = await import(
//         "@welldone-software/why-did-you-render"
//     );
//     console.log("Enabling whyDidYouRender...");
//     whyDidYouRender(React, {
//         trackAllPureComponents: true,
//         trackHooks: true,
//         logOnDifferentValues: true,
//         logOwnerReasons: true,
//         include: [/.*/],
//     });
// }
// const urlParams = new URLSearchParams(window.location.search);
// if (urlParams.get("wdyr")) {
//     loadWdyr();
// }
