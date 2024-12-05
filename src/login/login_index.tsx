// src/login/login_index.tsx
import { createRoot } from "react-dom/client";
import { CustomThemeProvider } from "@/ThemeProvider";
import Login from "./login";

const container = document.createElement("div");
document.body.appendChild(container);
const root = createRoot(container);

root.render(
    <CustomThemeProvider>
        <Login />
    </CustomThemeProvider>
);
