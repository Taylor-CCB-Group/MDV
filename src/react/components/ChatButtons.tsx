import { useCallback, useState } from "react";
import { useChartManager } from "../hooks";
import IconWithTooltip from "./IconWithTooltip";
import { ChatBubble as ChatBubbleIcon } from "@mui/icons-material";
import { CssBaseline, ThemeProvider, useTheme } from "@mui/material";
import useChat from "@/charts/dialogs/ChatAPI";
import ChatDialog from "@/charts/dialogs/ChatDialog";
import PopoutWindow from "./PopoutWindow";
import { ProjectProvider } from "@/modules/ProjectContext";

const ChatButtons = () => {
    const chatEnabled = useChartManager().config.chat_enabled;
    if (!chatEnabled) return null;

    const [open, setOpen] = useState(false);
    const [popout, setPopout] = useState(false);
    const theme = useTheme();

    // todo: come up with a different approach to only call it when the chat button is clicked
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

    const handlePopoutClose = useCallback(() => {
        setPopout(false);
    }, []);

    const handlePopoutOpen = useCallback(() => {
        setPopout(true);
    }, []);

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

    const onClose = useCallback(() => {
        setOpen(false);
    }, []);

    return (
        <>
            <IconWithTooltip tooltipText="Chat" onClick={() => setOpen(true)}>
                <ChatBubbleIcon />
            </IconWithTooltip>
            {!popout && <ChatDialog 
                open={open} 
                onClose={onClose} 
                onPopout={handlePopoutOpen} 
                {...chatDialogProps}
            />}
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
