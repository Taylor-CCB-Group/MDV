import { createRoot } from "react-dom/client";
   import Dashboard from "./Dashboard";
   import "./catalog_index.css";  // Make sure this line is present

   const container = document.createElement("div");
   document.body.appendChild(container);
   const root = createRoot(container);
   root.render(<Dashboard />);