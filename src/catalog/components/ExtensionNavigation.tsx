import { Launch } from "@mui/icons-material";
import { Box, Button } from "@mui/material";
import { buildApiUrl } from "@/utils/mdvRouting";
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
                <Button
                    variant={"outlined"}
                    key={item.id}
                    href={buildApiUrl(item.url)}
                    sx={{
                        borderWidth: 1.5,
                        borderRadius: 2,
                        px: 1.5,
                        py: 0.75,
                    }}
                    startIcon={<Launch />}
                >
                    {item.label}
                </Button>
            ))}
        </Box>
    );
}
