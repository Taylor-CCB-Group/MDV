import {
    Box,
    Button,
    Dialog,
    DialogContent,
    DialogTitle,
    Divider,
    Drawer,
    InputAdornment,
    List,
    ListItemButton,
    ListItemText,
    TextField,
    ThemeProvider,
    Typography,
    useTheme,
} from "@mui/material";
import type React from "react";
import Chatbot from "./ChatDialogComponent";
import { useChatLog } from "./ChatAPI";
import { Close, Launch, Search, ViewSidebar } from "@mui/icons-material";
import { useMemo, useState } from "react";
import IconWithTooltip from "@/react/components/IconWithTooltip";
import { usePopupWindow } from "./usePopupWindow";
import PopupChatDialog from "./PopupChatDialog";

export type ChatDialogProps = {
    open: boolean;
    setOpen: React.Dispatch<React.SetStateAction<boolean>>;
    fullscreen?: boolean;
};

const ChatDialog = ({ open, setOpen, fullscreen = false }: ChatDialogProps) => {
    const [drawerOpen, setDrawerOpen] = useState(true);
    const [search, setSearch] = useState("");
    // const [selectedChatId, setSelectedChatId] = useState();
    const handleClose = () => setOpen(false);
    const { chatLog } = useChatLog();
    // todo: Move it to a hook
    const filteredLog = useMemo(
        () => chatLog.filter((log) => log.query.toLowerCase().includes(search.toLowerCase())),
        [chatLog, search]
    );
    const theme = useTheme();

    const [isPopupMode, setIsPopupMode] = useState(false);
    const { openPopup, closePopup } = usePopupWindow({
      windowName: 'chat-dialog',
      onClose: () => {
        setIsPopupMode(false);
        setOpen(true);
      }
    });
    const handlePopout = () => {
      const popupContent = (
        <PopupChatDialog
          onClose={() => {
            closePopup();
            setIsPopupMode(false);
            setOpen(true);
          }}
        />
      );
      openPopup(popupContent, theme);
      setIsPopupMode(true);
      setOpen(false);
    };

    return (
        <Dialog open={open && !isPopupMode} onClose={handleClose} fullWidth maxWidth={"lg"} fullScreen={fullscreen}>
            <DialogTitle sx={{ px: 0, bgcolor: "var(--fade_background_color)" }}>
                <Box
                    sx={{
                        display: "flex",
                        alignItems: "center",
                    }}
                >
                    <IconWithTooltip
                        tooltipText={drawerOpen ? "Close Sidebar" : "Open Sidebar"}
                        onClick={() => setDrawerOpen(!drawerOpen)}
                        iconButtonProps={{
                            sx: {
                                marginLeft: 2,
                            },
                        }}
                    >
                        <ViewSidebar />
                    </IconWithTooltip>
                    <IconWithTooltip tooltipText={"Popout to a new window"} onClick={handlePopout}>
                        <Launch />
                    </IconWithTooltip>
                    <Typography variant="h6" sx={{ flexGrow: 1, textAlign: "center" }}>
                        {/* todo: Change this to the name of selected chat */}
                        Selected Chat
                    </Typography>

                    <IconWithTooltip
                        tooltipText={"Close"}
                        onClick={handleClose}
                        iconButtonProps={{
                            sx: {
                                marginRight: 2,
                            },
                        }}
                    >
                        <Close />
                    </IconWithTooltip>
                </Box>
            </DialogTitle>
            <DialogContent
                sx={{
                    height: "80vh",
                    padding: 0,
                    bgcolor: "var(--fade_background_color)",
                }}
                dividers
            >
                <Box sx={{ display: "flex", height: "100%" }}>
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
                                    value={search}
                                    onChange={(e) => setSearch(e.target.value)}
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
                                    {filteredLog.map((item) => (
                                        // todo: Add onclick to select chat id
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
                    <Box sx={{ flexGrow: 1, overflow: "hidden", pb: 2 }}>
                        <Chatbot />
                    </Box>
                </Box>
            </DialogContent>
        </Dialog>
    );
};

export default ChatDialog;
