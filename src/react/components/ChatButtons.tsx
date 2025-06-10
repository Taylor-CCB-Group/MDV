import { useCallback, useState } from "react";
import { useChartManager } from "../hooks";
import IconWithTooltip from "./IconWithTooltip";
import { ChatBubble as ChatBubbleIcon } from "@mui/icons-material";
import ChatButtonsComponent from "./ChatButtonsComponent";

const ChatButtons = () => {
    const chatEnabled = useChartManager().config.chat_enabled;
    const [open, setOpen] = useState(false);

    const onClose = useCallback(() => {
        setOpen(false);
    }, []);

    if (!chatEnabled) return null;

    return (
        <>
            <IconWithTooltip tooltipText="Chat" onClick={() => setOpen(true)}>
                <ChatBubbleIcon />
            </IconWithTooltip>
            {open && (
                <ChatButtonsComponent
                    open={open}
                    setOpen={setOpen}
                    onClose={onClose}
                />
            )}
        </>
    );
};

export default ChatButtons;
