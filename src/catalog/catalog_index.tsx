import React from "react";
import ReactDOM from "react-dom";
import App from "./catalog";

//console.log('this is the place to build the catalog/project browse interface...')
const root = document.createElement('div');
document.body.appendChild(root);
ReactDOM.render(React.createElement(App, null), root);