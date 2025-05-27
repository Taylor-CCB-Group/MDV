import {
    Box,
    Button,
    Dialog,
    DialogContent,
    DialogTitle,
    Divider,
    Drawer,
    IconButton,
    InputAdornment,
    List,
    ListItemButton,
    ListItemText,
    TextField,
    Tooltip,
} from "@mui/material";
import type React from "react";
import Chatbot from "./ChatDialogComponent";
import { useChatLog } from "./ChatAPI";
import { DialogCloseIconButton } from "@/catalog/ProjectRenameModal";
import { Search, ViewSidebar } from "@mui/icons-material";
import { useState } from "react";
import IconWithTooltip from "@/react/components/IconWithTooltip";

export type ChatDialogProps = {
    open: boolean;
    setOpen: React.Dispatch<React.SetStateAction<boolean>>;
};

const ChatDialog = ({ open, setOpen }: ChatDialogProps) => {
    const [drawerOpen, setDrawerOpen] = useState(true);
    const handleClose = () => setOpen(false);
    const { chatLog } = useChatLog();
    return (
        <Dialog open={open} onClose={handleClose} fullWidth maxWidth={"lg"}>
            <DialogTitle>
                <DialogCloseIconButton onClose={handleClose} />
            </DialogTitle>
            <DialogContent
                sx={{
                    height: "80vh",
                    padding: 0,
                }}
            >
                <Box sx={{ display: "flex", height: "100%" }}>
                    <IconWithTooltip
                        tooltipText={drawerOpen ? "Close Sidebar" : "Open Sidebar"}
                        onClick={() => setDrawerOpen(!drawerOpen)}
                        iconButtonProps={{
                            sx: {
                                position: "absolute",
                                top: 8,
                                left: 8,
                            },
                        }}
                    >
                        <ViewSidebar />
                    </IconWithTooltip>
                    {drawerOpen && (
                        <Drawer
                            open={drawerOpen}
                            variant="persistent"
                            anchor="left"
                            PaperProps={{
                                sx: {
                                    position: "relative",
                                    bgcolor: "unset",
                                },
                            }}
                            sx={{
                                width: 250,
                                height: "100%",
                            }}
                        >
                            <Box>
                                <Box sx={{ mt: 2, px: 1 }}>
                                    <Button variant="contained" fullWidth>
                                        New Chat
                                    </Button>
                                </Box>
                                <TextField
                                    fullWidth
                                    placeholder="Search History...."
                                    sx={{
                                        mt: 2,
                                        mb: 2,
                                        px: 1,
                                    }}
                                    slotProps={{
                                        input: {
                                            startAdornment: (
                                                <InputAdornment position="start">
                                                    <Search fontSize="small" />
                                                </InputAdornment>
                                            ),
                                        },
                                    }}
                                />
                                <Divider />
                                <List>
                                    {chatLog.map((item) => (
                                        <ListItemButton key={item.query}>
                                            <ListItemText
                                                primary={item.query}
                                                primaryTypographyProps={{
                                                    noWrap: true,
                                                    variant: "body1",
                                                    fontWeight: "bold",
                                                }}
                                            />
                                        </ListItemButton>
                                    ))}
                                </List>
                            </Box>
                        </Drawer>
                    )}
                    <Box sx={{ flexGrow: 1, overflow: "hidden", pt: 3, pb: 2 }}>
                        <Chatbot />
                    </Box>
                </Box>
            </DialogContent>
        </Dialog>
    );
};

export default ChatDialog;