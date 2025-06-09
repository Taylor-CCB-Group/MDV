import useChat from "@/charts/dialogs/ChatAPI";
import ChatDialog from "@/charts/dialogs/ChatDialog";
import { CssBaseline, ThemeProvider, useTheme } from "@mui/material";
import { useCallback, useState, type Dispatch, type SetStateAction } from "react";
import PopoutWindow from "./PopoutWindow";
import { ProjectProvider } from "@/modules/ProjectContext";

export type ChatButtonsComponentProps = {
    open: boolean,
    setOpen: Dispatch<SetStateAction<boolean>>;
};

const ChatButtonsComponent = ({open, setOpen}: ChatButtonsComponentProps) => {
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

    const onClose = useCallback(() => setOpen(false), [setOpen]);
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
    return (
        <>
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
    )
};

export default ChatButtonsComponent;