import { BotMessageSquare, SquareTerminal } from 'lucide-react';
import { MessageCircleQuestion, ThumbsUp, ThumbsDown, Star, NotebookPen, CircleAlert } from 'lucide-react';
import { type ChatProgress, type ChatMessage, navigateToView } from './ChatAPI';
import { useCallback, useEffect, useLayoutEffect, useMemo, useRef, useState } from 'react';
import JsonView from 'react18-json-view';
import ReactMarkdown from 'react-markdown';
import SyntaxHighlighter from 'react-syntax-highlighter';
import { dracula } from 'react-syntax-highlighter/dist/esm/styles/hljs';
import RobotPandaSVG from './PandaSVG';
import LinearProgress from '@mui/material/LinearProgress';
import { Box, Button, Divider, IconButton, InputAdornment, Skeleton, TextField } from '@mui/material';
import _ from 'lodash';
import { Check, ContentCopy, Clear as ClearIcon } from '@mui/icons-material';
import remarkGfm from 'remark-gfm';
import rehypeRaw from 'rehype-raw';
import rehypeSanitize from 'rehype-sanitize';


export type MessageType = {
    onClose: () => void;
    updateInput: (text: string) => void;
    suggestedQuestions: string[];
}

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

const SuggestedQuestions = ({ onSelect, suggestedQuestions }: { onSelect: (q: string) => void, suggestedQuestions: string[] }) => (
    <div className="flex flex-wrap gap-2 mt-4 mb-4 flex-col">
        <div className='font-bold'>Suggested Questions:</div>
        <div className='flex flex-wrap gap-2 mt-2'>{suggestedQuestions.map((q) => (
            <Button
                key={q}
                variant="outlined"
                onClick={() => onSelect(q)}
                size='large'
                sx={{ textTransform: 'none', borderRadius: 2, color: "inherit" }}
            >
                {q}
            </Button>
        ))}
        </div>
    </div>
);

const Message = ({ text: originalText, sender, view, onClose, error, updateInput, suggestedQuestions }: ChatMessage & MessageType) => {
    const isUser = sender === 'user';
    const [copied, setCopied] = useState(false);
    const pythonSections = extractPythonSections(originalText);
    
    let displayContent: any = originalText;
    try {
        displayContent = JSON.parse(originalText);
    } catch (e) {
        // Not JSON, keep as original text
    }

    const messageStyle = isUser ? 
        'bg-teal-200 self-end dark:bg-teal-900' :
        (error ?
            'bg-slate-200 dark:bg-slate-800 self-start border border-red-900':
            'bg-slate-200 dark:bg-slate-800 self-start'
        )
    
    const messageIcon = error ?
            <CircleAlert color='red' /> :
            isUser ? <MessageCircleQuestion className=''/> : <BotMessageSquare className='scale-x-[-1]' />;
    
    const handleCopy = async () => {
        try {
            const copyText = typeof displayContent === 'string' ? displayContent : JSON.stringify(displayContent, null, 2);
            await navigator.clipboard.writeText(copyText);
            setCopied(true);
            setTimeout(() => setCopied(false), 2000);
        } catch (err) {
            console.error("Failed to copy error details:", err);
        }
    };

    return (//setting `select-all` here doesn't help because * selector applies it to children, so we have custom class in tailwind theme
        <div className='selectable mt-4'>
            <div>{messageIcon}</div>
            <div className={`mb-2 p-4 rounded-lg ${messageStyle} relative markdown-body`}>
                {(sender === "bot" || error) && (
                    <IconButton
                        size="small"
                        aria-label="Copy button"
                        onClick={handleCopy}
                        sx={{ position: 'absolute', top: 2, right: 2 }}
                    >
                        {copied ? (
                            <Check fontSize='small' />
                        ) : (
                            <ContentCopy fontSize='small' />
                        )}
                    </IconButton>
                )}
                {typeof displayContent === 'object' && displayContent !== null ? (
                    <JsonView src={displayContent} />
                ) : (
                    <MessageMarkdown text={displayContent} />
                )}
            </div>
            {/* {pythonSections.map((section, index) => (
                <PythonCode key={index} code={section} />
            ))} */}

            {/* Show suggested questions only for the welcome message (sender: 'system') */}
            {sender === 'system' && suggestedQuestions.length > 0 && (
                <SuggestedQuestions suggestedQuestions={suggestedQuestions} onSelect={updateInput} />
            )}
            
            {/* Uncomment later and add logic for feedback buttons */}
            {/* {(sender === 'bot') && <MessageFeedback />} */}
            {view && 
                <Button 
                    variant="contained" 
                    color="primary" 
                    onClick={() => navigateToView(view, false, onClose)}
                >
                    Load view '{view}'...
                </Button>
            }
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

const MessageMarkdown = ({ text }: { text: string | any }) => {
    const markdown = typeof text === 'string' ? text : JSON.stringify(text, null, 2);
    // nb, I asked the bot for a markdown test, and a few things were in this markdown rendering
    // ~~strikethrough~~, tables, task-lists, and footnotes. Would be possible to use a plugin for that.
    // Also the example image, but that was probably because it was a bad link.
    return (
        <ReactMarkdown
            // biome-ignore lint/correctness/noChildrenProp: this is an issue with react-markdown, not our code
            children={markdown}
            remarkPlugins={[remarkGfm]}
            rehypePlugins={[rehypeRaw, rehypeSanitize]}
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
    const verbose = useMemo(() => 
        props.verboseProgress.map(s => s.substring(s.length-100)).join('\n'), 
    [props.verboseProgress]);

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

export type ChatBotProps = {
    messages: ChatMessage[];
    isSending: boolean;
    sendAPI: (input: string) => Promise<void>;
    requestProgress: ChatProgress | null;
    verboseProgress: string[];
    onClose: () => void;
    suggestedQuestions: string[];
};


const Chatbot = ({messages, isSending, sendAPI, requestProgress, verboseProgress, onClose, suggestedQuestions}: ChatBotProps) => {
    const [input, setInput] = useState<string>('');
    const messagesEndRef = useRef<HTMLDivElement>(null);
    // useCheckDataStore();

    const handleSend = useCallback(async () => {
        if (!input.trim()) return;
        setInput('');
        await sendAPI(input);
    }, [input, sendAPI]);

    const handleInputChange = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
        setInput(e.target.value);
    }, []);

    const updateInput = useCallback((text: string) => {
        setInput(text);
    }, []);

    const handleKeyPress = useCallback((e: React.KeyboardEvent<HTMLInputElement>) => {
        if (e.key === 'Enter') {
            handleSend();
        }
    }, [handleSend]);

    const handleClearInput = useCallback(() => {
        setInput("");
    }, []);

    const scrollToBottom = useCallback(() => {
        //! block: 'nearest' seems to solve issue with dialog header disappearing from top
        messagesEndRef.current?.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
    }, []);

    useLayoutEffect(() => {
        // Only scroll when there are new messages or progress updates
        if (messages || requestProgress) {
            scrollToBottom();
        }
    }, [messages, requestProgress, scrollToBottom]);
    
    return (
        <Box className="flex flex-col h-full mx-auto overflow-hidden">
            <Box className="flex-1 p-4 w-full overflow-y-auto" sx={{ bgcolor: "var(--background_color)"}}>
                {messages.map((message) => (
                    <Message key={`${message.id}-${message.sender}`} onClose={onClose} {...message} updateInput={updateInput} suggestedQuestions={suggestedQuestions} />
                ))}
                {isSending ?  
                    (requestProgress ? 
                        <Progress {...requestProgress} verboseProgress={verboseProgress} /> : 
                        <Skeleton variant='rectangular' sx={{ fontSize: "1.2rem", marginBottom: 2, marginTop: 2}} />
                    ) : (
                        <></>
                    )
                }
                {/* {
                isSending && 
                (<div className="animate-pulse flex justify-center p-4">{progressText}</div>)} */}
                <Box ref={messagesEndRef} />
            </Box>
             {/*<Box className='absolute opacity-10 pointer-events-none right-0'>
                <RobotPandaSVG />
            </Box> */}
            <Divider />
            <Box className="flex p-4 w-full">
                <TextField
                    type="text"
                    // disabled={isSending} //we can still type while it's processing
                    value={input}
                    onChange={handleInputChange}
                    onKeyDown={handleKeyPress}
                    placeholder="Type a message..."
                    fullWidth
                    sx={{
                        mr: 2,
                    }}
                    slotProps={{
                        input: {
                            endAdornment: (
                                <InputAdornment position="end">
                                    <IconButton onClick={handleClearInput}>
                                        <ClearIcon fontSize="small" />
                                    </IconButton>
                                </InputAdornment>
                            )
                        }
                    }}
                />
                <Button onClick={handleSend} disabled={isSending || !input} 
                className="p-2 bg-blue-500 text-white rounded-lg" variant='contained'>
                    Send
                </Button>
            </Box>
        </Box>
    );
};

export default Chatbot;