import { useProject } from "@/modules/ProjectContext";
import axios from "axios";
import { useEffect, useState } from "react";
import { z } from 'zod';

const responseSchema = z.object({
    message: z.string(),
    // viewName: z.string(),
});

// nb this is going to want to change a lot
// we want to know what version of code it based the training on as well as a bunch of other things
// and we probably want to be a bit flexible about what we accept and how to display it...
// Has the user said they like or dislike the response? Made any notes?
// we don't record that information yet - and it will mean a different way of interacting with the data
// not just logging new items, but updating existing ones - so message ID is important (currently just a random string in front-end)
const chatLogItemSchema = z.object({
    context: z.string(),
    query: z.string(),
    prompt_template: z.string(),
    response: z.string(),
});
const chatLogSchema = z.array(chatLogItemSchema);

export type ChatLogItem = z.infer<typeof chatLogItemSchema>;

export const useChatLog = () => {
    const { root } = useProject();
    const route = `${root}/chat_log.json`;
    const [chatLog, setChatLog] = useState<ChatLogItem[]>([]);
    const [isLoading, setIsLoading] = useState<boolean>(false);

    useEffect(() => {
        const fetchChatLog = async () => {
            setIsLoading(true);
            try {
                const response = await axios.get(route);
                const parsed = chatLogSchema.parse(response.data);
                setChatLog(parsed);
            } catch (error) {
                console.error('Error fetching chat log', error);
            }
            setIsLoading(false);
        };
        fetchChatLog();
        const interval = setInterval(fetchChatLog, 2000);
        return () => clearInterval(interval);
    }, []);

    return { chatLog, isLoading };
}

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