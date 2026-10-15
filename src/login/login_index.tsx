// src/login/login_index.tsx
import { createRoot } from "react-dom/client";
import { CustomThemeProvider } from "@/ThemeProvider";
import AccessDenied from "./AccessDenied";
import { getAccessDenial } from "./accessDenial";
import Login from "./login";

const container = document.createElement("div");
document.body.appendChild(container);
const root = createRoot(container);

// The server sends someone it turned away here with ?denied=<reason>, rather
// than to sign-in, which Auth0 would complete straight back into the denial.
const denial = getAccessDenial();

root.render(
    <CustomThemeProvider>
        {denial ? <AccessDenied denial={denial} /> : <Login />}
    </CustomThemeProvider>
);
