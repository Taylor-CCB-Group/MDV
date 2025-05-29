import { BotMessageSquare, SquareTerminal } from 'lucide-react';
import { MessageCircleQuestion, ThumbsUp, ThumbsDown, Star, NotebookPen } from 'lucide-react';
import useChat, { type ChatProgress, type ChatMessage, navigateToView } from './ChatAPI';
import { useCallback, useEffect, useRef, useState } from 'react';
import JsonView from 'react18-json-view';
import ReactMarkdown from 'react-markdown';
import SyntaxHighlighter from 'react-syntax-highlighter';
import { dracula } from 'react-syntax-highlighter/dist/esm/styles/hljs';
import RobotPandaSVG from './PandaSVG';
import LinearProgress from '@mui/material/LinearProgress';
import { Button } from '@mui/material';
import { useChartManager } from '@/react/hooks';
import { fetchJsonConfig } from '@/dataloaders/DataLoaderUtil';
import { useProject } from '@/modules/ProjectContext';
import type { DataSource } from '../charts';
import _ from 'lodash';
import ErrorDisplay from './DebugErrorComponent';


/**
 * prototype for new ways of handling project state...
 */
// function useCheckDataStore() {
//     const cm = useChartManager();
    
//     const { root } = useProject();
//     useEffect(() => {
//         // we want something that will respond via socket.io to changes in dataStore etc...
//         // this could be wrapped in a new more general-purpose hook, with an initial poll
//         // and then a socket.io connection to listen for changes
//         // but in either case we should be able to have a hook that returns the current state
//         // updated whenever it changes
//         const interval = setInterval(async () => {
//             // todo zod validation/schema
//             const config = await fetchJsonConfig(`${root}/datasources.json`, root) as DataSource[];
//             if (!config) {
//                 console.warn(`expected datasources.json, got ${config}`);
//                 return;
//             }
//             for (const ds of config) {
//                 const { name } = ds;
//                 // how ambitious do we want to be here in the short term?
//                 // we could just note that the state has changed, and therefore we need to refresh
//                 // but it would be much better to be able to update the state in place.
//                 // if we wanted to update size, that would be a problem, but we should be able to add
//                 // columns/links, new datasources... would be logical to have a DataStore.update method
//                 // but I'm inclined to have a function in a different module that is more strongly typed
//                 // and doesn't bloat the DataStore class
//                 if (!cm.dsIndex[name]) {
//                     console.log(`adding datasource ${name}`);
//                 }
//                 const old = cm.dsIndex[name]?.dataStore.config;
//                 if (!_.isEqual(old, ds)) {
//                     console.log(`updating datasource ${name}`);
//                 }
//             }
//         }, 1000);
//         return () => clearInterval(interval);
//     }, [cm, root]);
// }


const Message = ({ text, sender, view }: ChatMessage) => {
    const isUser = sender === 'user';
    const pythonSections = extractPythonSections(text);
    try {
        text = JSON.parse(text);
    } catch (e) {
    }
    const isError = sender === 'bot' && !view;
    return (//setting `select-all` here doesn't help because * selector applies it to children, so we have custom class in tailwind theme
        <div className='selectable'>
            {isUser ? <MessageCircleQuestion className=''/> : <BotMessageSquare className='scale-x-[-1]' />}
            <div className={`mb-2 p-4 rounded-lg ${
                isUser ? 'bg-teal-200 self-end dark:bg-teal-900' : 'bg-slate-200 dark:bg-slate-800 self-start'
                }`}>
                {/* <JsonView src={text} /> */}
                {isError ? <ErrorDisplay error={{ message: text }} /> : <MessageMarkdown text={text} />}
            </div>
            {/* {pythonSections.map((section, index) => (
                <PythonCode key={index} code={section} />
            ))} */}
            {(sender === 'bot') && <MessageFeedback />}
            {view && <Button variant="contained" color="primary" onClick={() => navigateToView(view)}>Load view '{view}'...</Button>}
        </div>
    );
}

const MessageFeedback = () => {
    const [isStarred, setIsStarred] = useState<boolean>(false);
    const [isLiked, setIsLiked] = useState<boolean>(false);
    const [isDisliked, setIsDisliked] = useState<boolean>(false);
    const [showNotes, setShowNotes] = useState<boolean>(false);
    const [notes, setNotes] = useState<string>('');
    return (
        <div className='flex justify-end'>
            <Star className='self-end scale-75' fill={isStarred ? 'var(--text_color)' : undefined}
            onClick={() => setIsStarred(!isStarred)}
            />
            <ThumbsUp className='self-end scale-75' fill={isLiked ? 'var(--text_color)' : undefined}
            onClick={() => {
                setIsLiked(!isLiked);
                setIsDisliked(false);
            }}
            />
            <ThumbsDown className='self-end scale-75' fill={isDisliked ? 'var(--text_color)' : undefined}
            onClick={() => {
                setIsDisliked(!isDisliked);
                setIsLiked(false);
            }}
            />
            <NotebookPen className='self-end scale-75' fill={showNotes ? 'var(--text_color)' : undefined} 
            onClick={() => setShowNotes(!showNotes)}
            />
            {showNotes && (
                <textarea
                    value={notes}
                    onChange={(e) => setNotes(e.target.value)}
                    placeholder='Add notes...'
                    className='p-2 m-2 border border-gray-300 rounded-lg'
                />
            )}
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
            {/* <pre className="text-sm font-mono text-gray-800 dark:text-gray-200">{code}</pre> */}
            <SyntaxHighlighter language="python" style={dracula}>
                {code}
            </SyntaxHighlighter>
        </div>
    );
}

const MessageMarkdown = ({ text }: { text: string }) => {
    const markdown = text;
    // nb, I asked the bot for a markdown test, and a few things were in this markdown rendering
    // ~~strikethrough~~, tables, task-lists, and footnotes. Would be possible to use a plugin for that.
    // Also the example image, but that was probably because it was a bad link.
    return (
        <ReactMarkdown
            // biome-ignore lint/correctness/noChildrenProp: this is an issue with react-markdown, not our code
            children={markdown}
            components={{
                code({ node, className, children, ...props }) {
                    const match = /language-(\w+)/.exec(className || '')
                    return match ? (
                        <>
                        <SquareTerminal onClick={() => alert(children)} /> {match[1]}:
                        <SyntaxHighlighter
                            className="rounded-lg border dark:border-gray-800 p-4 overflow-x-auto"
                            // biome-ignore lint/correctness/noChildrenProp: this is an issue with react-syntax-highlighter, not our code
                            children={String(children).replace(/\n$/, '')}
                            //@ts-ignore - not sure what's wrong here - maybe @types/react-syntax-highlighter needs updating?
                            style={dracula}
                            language={match[1]}
                            PreTag="div"
                            {...props}
                        />
                        </>
                    ) : (
                        <code className={className} {...props}>
                            {children}
                        </code>
                    )
                }
            }}
        />
    );
}

const Progress = (props: ChatProgress & {verboseProgress: string[]}) => {
    const verbose = props.verboseProgress.map(s => s.substring(s.length-100)).join('\n');
    if (props.progress >= 100) return null;
    return (
        <div className="p-4">
            <LinearProgress 
                variant="buffer"
                // consider some custom styling of the buffer
                value={props.progress}
                valueBuffer={props.progress + props.delta}
            />
            {props.message}
            <pre className="">{verbose}</pre>
        </div>
    );
}


const Chatbot = () => {
    const { messages, isSending, sendAPI, requestProgress, verboseProgress } = useChat();
    const [input, setInput] = useState<string>('');
    const messagesEndRef = useRef<HTMLDivElement>(null);
    // useCheckDataStore();

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
        //! block: 'nearest' seems to solve issue with dialog header disappearing from top
        messagesEndRef.current?.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
    }, []);

    useEffect(() => {
        messages;
        requestProgress;
        scrollToBottom();
    }, [messages, requestProgress, scrollToBottom]);
    
    return (
        <div className="flex flex-col h-full mx-auto overflow-hidden">
            <div className="flex-1 p-1 w-full overflow-y-auto">
                {messages.map((message) => (
                    <Message key={message.id} {...message} />
                ))}
                {requestProgress && <Progress {...requestProgress} verboseProgress={verboseProgress} />}
                {/* {
                isSending && 
                (<div className="animate-pulse flex justify-center p-4">{progressText}</div>)} */}
                <div ref={messagesEndRef} />
            </div>
            <div className='absolute opacity-10 pointer-events-none top-0 right-0'>
                <RobotPandaSVG />
            </div>
            <div className="flex p-4 border-t w-full border-gray-300">
                <input
                    type="text"
                    // disabled={isSending} //we can still type while it's processing
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