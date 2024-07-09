import { BotMessageSquare } from 'lucide-react';
import { MessageCircleQuestion } from 'lucide-react';
import useChat from './ChatAPI';
import { useCallback, useEffect, useRef, useState } from 'react';

const Message = ({ text, sender }: { text: string; sender: 'user' | 'bot' }) => {
    return (
        <div className={`mb-2 p-4 rounded-lg ${sender === 'user' ? 'bg-teal-200 self-end dark:bg-teal-900' : 'bg-slate-200 dark:bg-slate-800 self-start'}`}>
            {sender === 'user' ? <MessageCircleQuestion className=''/>
            : <BotMessageSquare className='scale-x-[-1]' />}
            {text}
        </div>
    );
}


const Chatbot = () => {
    const { messages, isSending, sendAPI } = useChat();
    const [input, setInput] = useState<string>('');
    const messagesEndRef = useRef<HTMLDivElement>(null);

    const handleSend = async () => {
        if (!input.trim()) return;
        await sendAPI(input);
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

    // biome-ignore lint/correctness/useExhaustiveDependencies: messages is needed for the effect to run
    useEffect(() => {
        scrollToBottom();
    }, [messages, scrollToBottom]);
    
    return (
        <div className="flex flex-col h-full mx-auto overflow-hidden">
            <div className="flex-1 p-4 overflow-y-auto">
                {messages.map((message) => (
                    <Message key={message.id} text={message.text} sender={message.sender} />
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