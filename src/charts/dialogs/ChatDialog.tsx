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
import Chatbot from "./ChatDialogComponent";
import type { ChatMessage, ChatProgress, ConversationMap } from "./ChatAPI";
import {
    Close as CloseIcon,
    Launch as LaunchIcon,
    Search as SearchIcon,
    ViewSidebar as ViewSidebarIcon,
} from "@mui/icons-material";
import { useMemo, useState, useRef, useCallback } from "react";
import IconWithTooltip from "@/react/components/IconWithTooltip";
import { useResizeDrawer } from "@/react/hooks";
import { Loader } from "@/react/components/ImportProjectDialog";

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
    suggestedQuestions: string[];
    onPopout?: () => void;
    isPopout?: boolean;
    fullscreen?: boolean;
    isLoading?: boolean;
    isLoadingInit?: boolean;
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
    suggestedQuestions,
    onPopout,
    isPopout,
    fullscreen = false,
    isLoading = false,
    isLoadingInit = false,
}: ChatDialogProps) => {
    const defaultDrawerWidth = 250;
    const minDrawerWidth = 180;
    const maxDrawerWidth = fullscreen ? 1000 : 600;
    const [drawerOpen, setDrawerOpen] = useState(true);
    const [search, setSearch] = useState("")
    const dialogRef = useRef<HTMLDivElement>(null);
    const { 
        drawerWidth,
        onMouseDown,
    } = useResizeDrawer(
        dialogRef, 
        defaultDrawerWidth, 
        minDrawerWidth, 
        maxDrawerWidth
    );

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
            ref={dialogRef}
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
                        <>
                        <Drawer
                            open={drawerOpen}
                            variant="persistent"
                            anchor="left"
                            PaperProps={{
                                sx: {
                                    position: "relative",
                                    bgcolor: "unset",
                                    width: drawerWidth,
                                    minWidth: minDrawerWidth,
                                    maxWidth: maxDrawerWidth,
                                },
                            }}
                            sx={{
                                width: drawerWidth,
                                height: "100%",
                            }}
                        >
                            <Box sx={{ position: 'relative', height: '100%' }}>
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
                                    {isLoading || isLoadingInit ? (
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
                        <Box
                            sx={{
                                // position: 'absolute',
                                top: 0,
                                left: 0,
                                width: 8,
                                height: '100%',
                                cursor: 'ew-resize',
                                zIndex: 10,
                                backgroundColor: 'rgba(0,0,0,0.05)',
                                '&:hover': { background: 'rgba(0,0,0,0.15)' },
                            }}
                            onMouseDown={onMouseDown}
                        />
                        </>
                    )}
                    <Box sx={{ flexGrow: 1, overflow: "hidden", pb: 2 }}>
                        {isLoadingInit ? 
                            (<Loader />) : 
                            (<Chatbot
                                messages={messages}
                                isSending={isSending}
                                requestProgress={requestProgress}
                                sendAPI={sendAPI}
                                verboseProgress={verboseProgress}
                                onClose={onClose}
                                suggestedQuestions={suggestedQuestions}
                            />)
                            }
                    </Box>
                </Box>
            </DialogContent>
        </Dialog>
    );
};

export default ChatDialog;
