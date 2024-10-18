import { createRoot } from "react-dom/client";
import Dashboard from "./Dashboard";
import "./catalog_index.css";
import { CustomThemeProvider } from "@/ThemeProvider";

const container = document.createElement("div");
document.body.appendChild(container);
const root = createRoot(container);

root.render(
  <CustomThemeProvider>
    <Dashboard />
  </CustomThemeProvider>
);