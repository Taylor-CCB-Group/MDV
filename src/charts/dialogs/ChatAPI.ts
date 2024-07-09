import { useProject } from "@/modules/ProjectContext";
import axios from "axios";
import { useState } from "react";
import { z } from 'zod';

const responseSchema = z.object({
    message: z.string(),
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



const sendMessage = async (message: string, route = '/chat') => {
    const response = await axios.post<ChatResponse>(route, { message });
    const parsed = responseSchema.parse(response.data); // may throw an error if the response is not valid
    return parsed;
};


const useChat = () => {
    const { root } = useProject();
    const route = `${root}/chat`;
    const [messages, setMessages] = useState<Message[]>([]);
    const [isSending, setIsSending] = useState<boolean>(false);

    const appendMessage = (message: string, sender: 'bot' | 'user') => {
        const msg = { text: message, sender, id: generateId() };
        setMessages((prevMessages) => [...prevMessages, msg]);
    };
    
    const sendAPI = async (input: string) => {
        if (!input.trim()) return;

        const userMessage: Message = { text: input, sender: 'user', id: generateId() };
        appendMessage(input, 'user');

        try {
            setIsSending(true);
            const response = await sendMessage(input, route);
            appendMessage(response.message, 'bot');
        } catch (error) {
            appendMessage(`Error: ${error}`, 'bot');
        }
        setIsSending(false);
        setIsSending(false);
    };
    return { messages, isSending, sendAPI };
}

export default useChat;