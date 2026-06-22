import { createRoot } from "react-dom/client";
import { CustomThemeProvider } from "@/ThemeProvider";
import AdminApp from "./AdminApp";
import "./admin.css";

const container =
    document.getElementById("admin-root") ?? document.body.appendChild(document.createElement("div"));

createRoot(container).render(
    <CustomThemeProvider>
        <AdminApp />
    </CustomThemeProvider>,
);
