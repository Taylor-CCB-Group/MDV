import { useCallback, useState } from "react";
import { useChartManager } from "../hooks";
import IconWithTooltip from "./IconWithTooltip";
import { ChatBubble as ChatBubbleIcon } from "@mui/icons-material";
import { CssBaseline, ThemeProvider, useTheme } from "@mui/material";
import useChat from "@/charts/dialogs/ChatAPI";
import ChatDialog from "@/charts/dialogs/ChatDialog";
import PopoutWindow from "./PopoutWindow";
import { ProjectProvider } from "@/modules/ProjectContext";
import type React from 'react';

interface ChatProviderProps {
    open: boolean;
    popout: boolean;
    onClose: () => void;
    onPopout: () => void;
    onPopoutClose: () => void;
    theme: any;
}

const ChatProvider: React.FC<ChatProviderProps> = ({
    open,
    popout,
    onClose,
    onPopout,
    onPopoutClose,
    theme
}) => {
    // useChat is only called when this component is mounted (after first button click)
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
        isLoadingInit,
        suggestedQuestions,
    } = useChat();

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
        isLoadingInit,
        suggestedQuestions,
    };

    const dialog = (
        <ChatDialog
            open={popout ? true : open}
            onClose={popout ? onPopoutClose : onClose}
            onPopout={!popout ? onPopout : undefined}
            fullscreen={popout}
            isPopout={popout}
            {...chatDialogProps}
        />
    );

    return (
        <>
            {popout ? (
                <PopoutWindow onClose={onPopoutClose}>
                    <ProjectProvider>
                        <ThemeProvider theme={theme}>
                            <CssBaseline />
                            {dialog}
                        </ThemeProvider>
                    </ProjectProvider>
                </PopoutWindow>
                ) : dialog
            }
        </>
    );
};

const ChatButtons = () => {
    //! this property doesn't change, so it's ok to return early which would otherwise violate the rule of hooks
    const chatEnabled = useChartManager().config.chat_enabled;
    if (!chatEnabled) return null;

    const [open, setOpen] = useState(false);
    const [popout, setPopout] = useState(false);
    const [chatInitialized, setChatInitialized] = useState(false);
    const theme = useTheme();

    const handlePopoutClose = useCallback(() => {
        setPopout(false);
    }, []);

    const handlePopoutOpen = useCallback(() => {
        setPopout(true);
    }, []);

    // Chat initialization happens on first button click - this is when useChat() gets called
    const handleChatOpen = useCallback(() => {
        if (!chatInitialized) {
            setChatInitialized(true);
        }
        setOpen(true);
    }, [chatInitialized]);

    const onClose = useCallback(() => {
        setOpen(false);
    }, []);

    return (
        <>
            <IconWithTooltip tooltipText="Chat" onClick={handleChatOpen}>
                <ChatBubbleIcon />
            </IconWithTooltip>

            {/* ChatProvider (and useChat) only mounts after first button click */}
            {chatInitialized && (
                <ChatProvider
                    open={open}
                    popout={popout}
                    onClose={onClose}
                    onPopout={handlePopoutOpen}
                    onPopoutClose={handlePopoutClose}
                    theme={theme}
                />
            )}
        </>
    );
};

export default ChatButtons;