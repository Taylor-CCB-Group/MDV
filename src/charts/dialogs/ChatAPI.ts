import { useProject } from "@/modules/ProjectContext";
import axios from "axios";
import { useEffect, useState } from "react";
import { z } from 'zod';

const responseSchema = z.object({
    message: z.string(),
    // viewName: z.string(),
});

type ChatResponse = z.infer<typeof responseSchema>;

type Message = {
    text: string;
    sender: 'user' | 'bot';
    id: string;
};

function generateId() {
    return Math.random().toString(36).substring(7);
}
/** viewName could be a prop of ProjectProvider, but currently not cleanly reactive */
function getViewName(): string | null {
    const urlParams = new URLSearchParams(window.location.search);
    return urlParams.get('view');
}


const sendMessage = async (message: string, route = '/chat') => {
    // we should send information about the context - in particular, which view we're in
    const response = await axios.post<ChatResponse>(route, { message });
    const parsed = responseSchema.parse(response.data); // may throw an error if the response is not valid
    return parsed;
};

const DefaultMessage: Message = {
    text: 'Hello! How can I help you?',
    sender: 'bot',
    id: generateId(),
}

const useChat = () => {
    const { root } = useProject(); //todo add viewName to ProjectProvider
    const route = `${root}/chat`;
    const routeInit = `${root}/chat_init`;
    const [messages, setMessages] = useState<Message[]>([]);
    const [isSending, setIsSending] = useState<boolean>(false);
    const [isInit, setIsInit] = useState<boolean>(false);

    useEffect(() => {
        const chatInit = async () => {
            setIsSending(true);
            try {
                const response = await sendMessage('', routeInit);
                setMessages([{ text: response.message, sender: 'bot', id: generateId() }]);
            } catch (error) {
                console.error('Error sending welcome message', error);
            }
            setIsSending(false);
            setIsInit(true);
        };
        if (!isSending && !isInit) chatInit();
    }, [isSending]);

    const appendMessage = (message: string, sender: 'bot' | 'user') => {
        const msg = { text: message, sender, id: generateId() };
        setMessages((prevMessages) => [...prevMessages, msg]);
    };
    
    const sendAPI = async (input: string) => {
        if (!input.trim()) return;
        const viewName = getViewName();
        console.log(`sending chat '${input}' from '${viewName}'`)
        appendMessage(input, 'user');

        try {
            setIsSending(true);
            const response = await sendMessage(input, route);
            appendMessage(response.message, 'bot');
        } catch (error) {
            appendMessage(`Error: ${error}`, 'bot');
        }
        setIsSending(false);
    };
    return { messages, isSending, sendAPI };
}

export default useChat;