import { createRoot } from "react-dom/client";
import { CustomThemeProvider } from "@/ThemeProvider";
import AdminApp from "./AdminApp";
import "./admin.css";

const container = document.getElementById("admin-root") ?? document.createElement("div");
container.id = "admin-root";
if (!container.parentElement) {
    document.body.appendChild(container);
}

createRoot(container).render(
    <CustomThemeProvider>
        <AdminApp />
    </CustomThemeProvider>,
);
