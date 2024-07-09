import { BotMessageSquare } from 'lucide-react';
import { MessageCircleQuestion } from 'lucide-react';


import { sendMessage } from './ChatAPI';
import { useCallback, useEffect, useRef, useState } from 'react';

type Message = {
    text: string;
    sender: 'user' | 'bot';
    id: string;
};

function generateId() {
    return Math.random().toString(36).substring(7);
}



const Chatbot = () => {
    const [messages, setMessages] = useState<Message[]>([]);
    const [input, setInput] = useState<string>('');
    const [isSending, setIsSending] = useState<boolean>(false);
    const messagesEndRef = useRef<HTMLDivElement>(null);

    const handleSend = async () => {
        if (!input.trim()) return;

        const userMessage: Message = { text: input, sender: 'user', id: generateId() };
        setMessages((prevMessages) => [...prevMessages, userMessage]);

        try {
            setIsSending(true);
            const response = await sendMessage(input);
            const botMessage: Message = { text: response.message, sender: 'bot', id: generateId() };
            setMessages((prevMessages) => [...prevMessages, botMessage]);
        } catch (error) {
            const errorMessage: Message = { text: `Error: ${error}`, sender: 'bot', id: generateId() };
            setMessages((prevMessages) => [...prevMessages, errorMessage]);
        }
        setIsSending(false);

        setInput('');
    };

    const handleInputChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        setInput(e.target.value);
    };

    const handleKeyPress = (e: React.KeyboardEvent<HTMLInputElement>) => {
        if (e.key === 'Enter') {
            handleSend();
        }
    };
    const scrollToBottom = useCallback(() => {
        messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
    }, []);

    // biome-ignore lint/correctness/useExhaustiveDependencies: messages is 
    useEffect(() => {
        scrollToBottom();
    }, [messages, scrollToBottom]);
    
    return (
        <div className="flex flex-col h-full mx-auto overflow-hidden">
            <div className="flex-1 p-4 overflow-y-auto">
                {messages.map((message) => (
                    <div key={message.id} className={`mb-2 p-2 rounded-lg ${message.sender === 'user' ? 'bg-teal-200 self-end dark:bg-teal-900' : 'bg-slate-200 dark:bg-slate-800 self-start'}`}>
                        {message.sender === 'user' ? <MessageCircleQuestion /> : <BotMessageSquare />}
                        {message.text}
                    </div>
                ))}
                <div ref={messagesEndRef} />
            </div>
            <div className="flex p-4 border-t border-gray-300">
                <input
                    type="text"
                    disabled={isSending}
                    value={input}
                    onChange={handleInputChange}
                    onKeyDown={handleKeyPress}
                    placeholder="Type a message..."
                    className="flex-1 p-2 border border-gray-300 rounded-lg mr-2"
                />
                <button type="submit" onClick={handleSend} disabled={isSending} 
                className="p-2 bg-blue-500 text-white rounded-lg">
                    Send
                </button>
            </div>
        </div>
    );
};

export default Chatbot;