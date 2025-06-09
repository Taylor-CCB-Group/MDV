import { useCallback, useEffect, useState } from "react";
import { useChartManager } from "../hooks";
import { CssBaseline, ThemeProvider, useTheme } from "@mui/material";
import useChat from "@/charts/dialogs/ChatAPI";
import IconWithTooltip from "./IconWithTooltip";
import { ChatBubble as ChatBubbleIcon } from "@mui/icons-material";
import ChatDialog from "@/charts/dialogs/ChatDialog";
import { ProjectProvider } from "@/modules/ProjectContext";
import PopoutWindow from "./PopoutWindow";

const ChatButtons = () => {
    const chatEnabled = useChartManager().config.chat_enabled;
    const [open, setOpen] = useState(false);
    const [popout, setPopout] = useState(false);
    const theme = useTheme();

    const {
        messages,
        isSending,
        sendAPI,
        requestProgress,
        verboseProgress,
        isChatLogLoading,
        startNewConversation,
        switchConversation,
        conversationMap,
        conversationId,
    } = useChat();

    const onClose = useCallback(() => setOpen(false), []);
    const handlePopoutClose = useCallback(() => {
        setPopout(false);
    }, []);

    const handlePopoutOpen = useCallback(() => {
        setPopout(true);
    }, []);

    if (!chatEnabled) return null;

    const chatDialogProps = {
        messages,
        isSending,
        requestProgress,
        sendAPI,
        verboseProgress,
        isLoading: isChatLogLoading,
        startNewConversation,
        switchConversation,
        conversationMap,
        conversationId,
    };

    return (
        <>
            <IconWithTooltip tooltipText="Chat" onClick={() => setOpen(true)}>
                <ChatBubbleIcon />
            </IconWithTooltip>
            {!popout && <ChatDialog open={open} onClose={onClose} onPopout={handlePopoutOpen} {...chatDialogProps} />}
            {popout && (
                <PopoutWindow onClose={handlePopoutClose}>
                    <ProjectProvider>
                        <ThemeProvider theme={theme}>
                            <CssBaseline />
                            <ChatDialog
                                open={true}
                                onClose={handlePopoutClose}
                                isPopout
                                fullscreen
                                {...chatDialogProps}
                            />
                        </ThemeProvider>
                    </ProjectProvider>
                </PopoutWindow>
            )}
        </>
    );
};

export default ChatButtons;
