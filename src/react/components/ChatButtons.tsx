import { useCallback, useState } from "react";
import { useChartManager } from "../hooks";
import { CssBaseline, ThemeProvider, useTheme } from "@mui/material";
import useChat, { useChatLog } from "@/charts/dialogs/ChatAPI";
import ChatLogDialog from "@/charts/dialogs/ChatLogDialog";
import IconWithTooltip from "./IconWithTooltip";
import { ChatBubble as ChatBubbleIcon } from "@mui/icons-material";
import ChatLogIcon from "@mui/icons-material/WebStories";
import ChatDialog from "@/charts/dialogs/ChatDialog";
import { ProjectProvider } from "@/modules/ProjectContext";
import PopoutWindow from "./PopoutWindow";

const ChatButtons = () => {
    // very basic check to see if chat is enabled
    const chatEnabled = useChartManager().config.chat_enabled;
    const [open, setOpen] = useState(false);
    const [popout, setPopout] = useState(false);
    const theme = useTheme();

    const { chatLog } = useChatLog();
    const { messages, isSending, sendAPI, requestProgress, verboseProgress } = useChat();

    const onClose = () => setOpen(false);
    const handlePopoutClose = useCallback(() => {
        setPopout(false);
    }, []);

    const handlePopoutOpen = useCallback(() => {
        setPopout(true);
    }, []);

    const handleChatLogButtonClick = () => {
        new ChatLogDialog();
    };
    if (!chatEnabled) return null;
    return (
        <>
            <IconWithTooltip tooltipText="Chat" onClick={() => setOpen(true)}>
                <ChatBubbleIcon />
            </IconWithTooltip>
            <IconWithTooltip tooltipText="Chat Log" onClick={handleChatLogButtonClick}>
                <ChatLogIcon />
            </IconWithTooltip>
            {!popout && (
                <ChatDialog
                    open={open}
                    onClose={onClose}
                    onPopout={handlePopoutOpen}
                    chatLog={chatLog}
                    messages={messages}
                    isSending={isSending}
                    requestProgress={requestProgress}
                    sendAPI={sendAPI}
                    verboseProgress={verboseProgress}
                />
            )}
            {popout && (
                <PopoutWindow onClose={handlePopoutClose}>
                    <ProjectProvider>
                        <ThemeProvider theme={theme}>
                            <CssBaseline />
                            <ChatDialog
                                open={true}
                                onClose={handlePopoutClose}
                                chatLog={chatLog}
                                messages={messages}
                                isSending={isSending}
                                requestProgress={requestProgress}
                                sendAPI={sendAPI}
                                verboseProgress={verboseProgress}
                                isPopout
                                fullscreen
                            />
                        </ThemeProvider>
                    </ProjectProvider>
                </PopoutWindow>
            )}
        </>
    );
};

export default ChatButtons;