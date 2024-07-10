import { BotMessageSquare, SquareTerminal } from 'lucide-react';
import { MessageCircleQuestion } from 'lucide-react';
import useChat from './ChatAPI';
import { useCallback, useEffect, useRef, useState } from 'react';
import JsonView from 'react18-json-view';
import Markdown from 'react-markdown';

const Message = ({ text, sender }: { text: string; sender: 'user' | 'bot' }) => {
    const isUser = sender === 'user';
    const pythonSections = extractPythonSections(text);
    try {
        text = JSON.parse(text);
    } catch (e) {
    }
    return (
        <div>
            {isUser ? <MessageCircleQuestion className=''/> : <BotMessageSquare className='scale-x-[-1]' />}
            <div className={`mb-2 p-4 rounded-lg ${
                isUser ? 'bg-teal-200 self-end dark:bg-teal-900' : 'bg-slate-200 dark:bg-slate-800 self-start'
                }`}>
                {/* <JsonView src={text} /> */}
                <MessageMarkdown text={text} />
            </div>
            {pythonSections.map((section, index) => (
                <PythonCode key={index} code={section} />
            ))}
        </div>
    );
}
const MessageJson = ({ text }: { text: string }) => {
    try {
        text = JSON.parse(text);
    } catch (e) {
    }
    return (
        <JsonView src={text} />
    );
}
function extractPythonSections(responseText: string) {
    // Regular expression to match Python code blocks 
    const pythonCodeRegex = /```python([\s\S]*?)```/g;
    let matches; const pythonSections = [];
    // Find all matches 
    while ((matches = pythonCodeRegex.exec(responseText)) !== null) {
        // Extract the code block without the backticks and "python" keyword 
        pythonSections.push(matches[1].trim());
    }
    return pythonSections;
}

const PythonCode = ({ code }: { code: string }) => {
    return (
        <div className="p-4 bg-gray-200 dark:bg-gray-800 w-fit mb-4 rounded-lg">
            <SquareTerminal />
            <pre className="text-sm font-mono text-gray-800 dark:text-gray-200">{code}</pre>
        </div>
    );
}

const MessageMarkdown = ({ text }: { text: string }) => {
    // nb, I asked the bot for a markdown test, and a few things were in this markdown rendering
    // ~~strikethrough~~, tables, task-lists, and footnotes. 
    // Also the example image, but that was probably because it was a bad link.
    return (
        <Markdown>{text}</Markdown>
    );
}

const Chatbot = () => {
    const { messages, isSending, sendAPI } = useChat();
    const [input, setInput] = useState<string>('');
    const messagesEndRef = useRef<HTMLDivElement>(null);

    const handleSend = async () => {
        if (!input.trim()) return;
        setInput('');
        await sendAPI(input);
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
            <div className="flex-1 p-1 overflow-y-auto">
                {messages.map((message) => (
                    <Message key={message.id} text={message.text} sender={message.sender} />
                ))}
                <div ref={messagesEndRef} />
            </div>
            <div className="flex p-4 border-t w-full border-gray-300">
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