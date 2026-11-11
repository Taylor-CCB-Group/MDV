import { createRoot } from "react-dom/client";
import Dashboard from "./Dashboard";
import "./catalog_index.css";
import { CustomThemeProvider } from "@/ThemeProvider";
import { offerOwnTabIfFramed } from "@/utils/ownTabNotice";
import { PermissionsProvider } from "./PermissionsContext";

offerOwnTabIfFramed();

const container = document.createElement("div");
document.body.appendChild(container);
const root = createRoot(container);

root.render(
    <CustomThemeProvider>
        <PermissionsProvider>
            <Dashboard />
        </PermissionsProvider>
    </CustomThemeProvider>,
);
