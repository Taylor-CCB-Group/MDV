import axios from "axios";
import { useState } from "react";
import { z } from 'zod';

const responseSchema = z.object({
    message: z.string(),
});

type ChatResponse = z.infer<typeof responseSchema>;

export const sendMessage = async (message: string) => {
    const response = await axios.post<ChatResponse>('/chat', { message });
    const parsed = responseSchema.parse(response.data); // may throw an error if the response is not valid
    return parsed;
};


// const useChat = () => {
//     const [messages, setMessages] = useState<string[]>([]);
    
//     const appendMessage = (message: string) => {
//         setMessages((prevMessages) => [...prevMessages, message]);
//     };
    
//     return { messages, appendMessage };
// }

// export default useChat;