import { buildApiUrl } from "@/utils/mdvRouting";
import { Box, Button } from "@mui/material";
import useExtensionNavigation from "../hooks/useExtensionNavigation";
import { Launch } from "@mui/icons-material";


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
                    color="inherit" 
                    sx={{ 
                        borderWidth: 1,
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
