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
    Typography,
} from "@mui/material";
import Chatbot, { type ChatBotProps } from "./ChatDialogComponent";
import type { ChatLogItem, ChatMessage, ChatProgress, ConversationLog, ConversationMap } from "./ChatAPI";
import {
    Close as CloseIcon,
    Launch as LaunchIcon,
    Search as SearchIcon,
    ViewSidebar as ViewSidebarIcon,
} from "@mui/icons-material";
import { useEffect, useMemo, useState } from "react";
import IconWithTooltip from "@/react/components/IconWithTooltip";

export type ChatDialogProps = {
    open: boolean;
    onClose: () => void;
    messages: ChatMessage[];
    isSending: boolean;
    sendAPI: (input: string) => Promise<void>;
    requestProgress: ChatProgress | null;
    verboseProgress: string[];
    startNewConversation: () => Promise<void>;
    switchConversation: (id: string) => void;
    conversationMap: ConversationMap;
    conversationId: string;
    onPopout?: () => void;
    isPopout?: boolean;
    fullscreen?: boolean;
    isLoading?: boolean;
};

const ChatDialog = ({
    open,
    onClose,
    messages,
    isSending,
    sendAPI,
    requestProgress,
    verboseProgress,
    startNewConversation,
    switchConversation,
    conversationMap,
    conversationId,
    onPopout,
    isPopout,
    fullscreen = false,
    isLoading = false,
}: ChatDialogProps) => {
    const [drawerOpen, setDrawerOpen] = useState(true);
    const [search, setSearch] = useState("");

    const filteredLog = useMemo(() => {
        const filteredConversations: ConversationMap = {};
        Object.entries(conversationMap)
            .slice()
            .reverse()
            .forEach(([convId, conv]) => {
                const found = conv.logText.toLowerCase().includes(search.toLowerCase());
                if (found) {
                    filteredConversations[convId] = conv;
                }
            });
        return filteredConversations;
    }, [conversationMap, search]);

    return (
        <Dialog
            open={open}
            onClose={onClose}
            fullWidth
            maxWidth={"lg"}
            fullScreen={fullscreen}
            disablePortal={isPopout}
        >
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
                        <ViewSidebarIcon />
                    </IconWithTooltip>
                    {!isPopout && (
                        <IconWithTooltip tooltipText={"Popout to a new window"} onClick={() => onPopout?.()}>
                            <LaunchIcon />
                        </IconWithTooltip>
                    )}
                    <Typography variant="h6" sx={{ flexGrow: 1, textAlign: "center" }}>
                        Chat MDV
                    </Typography>

                    <IconWithTooltip
                        tooltipText={"Close"}
                        onClick={onClose}
                        iconButtonProps={{
                            sx: {
                                marginRight: 2,
                            },
                        }}
                    >
                        <CloseIcon />
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
                                    <Button variant="contained" fullWidth onClick={startNewConversation}>
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
                                                    <SearchIcon fontSize="small" />
                                                </InputAdornment>
                                            ),
                                        },
                                    }}
                                />
                                <Divider />
                                <List>
                                    {isLoading ? (
                                        <ListItemButton>
                                            <ListItemText primary="Loading chat history..." />
                                        </ListItemButton>
                                    ) : (
                                        Object.entries(filteredLog).map(([convId, conv], index) => (
                                            <ListItemButton
                                                key={convId}
                                                onClick={() => switchConversation(convId)}
                                                selected={conversationId === convId}
                                            >
                                                <ListItemText
                                                    primary={conv.logText}
                                                    secondary={`${conv.logLength} messages`}
                                                    primaryTypographyProps={{
                                                        noWrap: true,
                                                        variant: "body1",
                                                        fontWeight: "bold",
                                                    }}
                                                />
                                            </ListItemButton>
                                        ))
                                    )}
                                </List>
                            </Box>
                        </Drawer>
                    )}
                    <Box sx={{ flexGrow: 1, overflow: "hidden", pb: 2 }}>
                        <Chatbot
                            messages={messages}
                            isSending={isSending}
                            requestProgress={requestProgress}
                            sendAPI={sendAPI}
                            verboseProgress={verboseProgress}
                        />
                    </Box>
                </Box>
            </DialogContent>
        </Dialog>
    );
};

export default ChatDialog;
