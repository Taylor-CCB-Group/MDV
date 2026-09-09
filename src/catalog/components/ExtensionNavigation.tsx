import { buildApiUrl } from "@/utils/mdvRouting";
import { Box, Button } from "@mui/material";
import useExtensionNavigation from "../hooks/useExtensionNavigation";

export default function ExtensionNavigation() {
    const items = useExtensionNavigation();

    if (items.length === 0) {
        return null;
    }

    return (
        <Box
            component="nav"
            aria-label="Enabled extensions"
            sx={{ display: "flex", alignItems: "center", gap: 1, ml: 2 }}
        >
            {items.map((item) => (
                <Button key={item.id} href={buildApiUrl(item.url)} color="inherit" sx={{ textTransform: "none" }}>
                    {item.label}
                </Button>
            ))}
        </Box>
    );
}
