import { useProject } from "@/modules/ProjectContext";
import axios from "axios";
import { useCallback, useEffect, useMemo, useState } from "react";
import * as z from 'zod/v4';
import { useQuery } from '@tanstack/react-query';
import { useChartManager, useViewManager } from "@/react/hooks";

const completedChatResponseSchema = z.object({
    message: z.string(),
    /** 
     * This should be a string that identifies newly created views.
     */
    view: z.string().optional().nullable().describe('The name of the newly created view, if a view was created'),
    // timestamp: z.string(),
    error: z.boolean().optional().nullable(),
});
const chatInitResponseSchema = z.object({
    message: z.string(),
    suggested_questions: z.array(z.string()).optional().nullable(),
    error: z.boolean().optional().nullable(),
});
type ChatInitResponse = z.infer<typeof chatInitResponseSchema>;

const chatProgressSchema = z.object({
    id: z.string().describe('The ID of the message this progress is associated with'),
    progress: z.number(),
    delta: z.number(),
    message: z.string().optional().nullable(),
});

export type ChatProgress = z.infer<typeof chatProgressSchema>;

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
    // id: z.string(),
    //nullable as well as optional because of some data that had null during development...
    //the offending code was never merged so that is just for testing
    //these fields are expected to always be present on actual 
    conversation_id: z.string().optional().nullable(),
    timestamp: z.string().optional().nullable(),
    view_name: z.string().optional().nullable(),
    error: z.boolean().optional().nullable(),
});
// this is possible - may consider it at some point.
// .transform((data) => ({
//     convesationId: data.conversation_id,
// }));
const chatLogSchema = z.array(chatLogItemSchema);
// z.toJSONSchema(chatLogSchema);

export type ChatLogItem = z.infer<typeof chatLogItemSchema>;

type ChatResponse = z.infer<typeof completedChatResponseSchema>;

export type ChatMessage = {
    text: string;
    view?: string;
    sender: 'user' | 'bot' | 'system';
    id: string;
    conversationId: string;
    error?: boolean;
};


export type ConversationLog = {
    logText: string,
    logLength: number,
    messages: ChatMessage[],
}

export type ConversationMap = Record<string, ConversationLog>;

function generateId() {
    return Math.random().toString(36).substring(7);
}

function generateConversationId() {
    return `conv-${Date.now()}-${Math.random().toString(36).substring(7)}`;
}

/** viewName could be a prop of ProjectProvider, but currently not cleanly reactive */
function getViewName(): string | null {
    const urlParams = new URLSearchParams(window.location.search);
    return urlParams.get('view');
}

// this is now only used for the initial message, meaning that we can change the response format
// to include suggested questions.
const sendChatInitHttp = async (message: string, id: string, route: string, conversationId: string) => {
    // we should send information about the context - in particular, which view we're in
    // could consider a streaming response here rather than socket
    const response = await axios.post<ChatInitResponse>(route, { message, id, conversation_id: conversationId });
    const parsed = chatInitResponseSchema.parse(response.data); // may throw an error if the response is not valid
    return parsed;
};

const sendMessageSocket = async (message: string, id: string, _routeUnused: string, conversationId: string, time?: number) => {
    // consider refactoring to be a generator function, so that we can yield progress updates
    const socket = window.mdv.chartManager.ipc?.socket;
    if (!socket) return;

    socket.emit('chat_request', { message, id, conversation_id: conversationId });
    const response = await new Promise<ChatResponse>((resolve, reject) => {
        const timeout = time ? setTimeout(() => {
            socket.off('chat_response', onChatResponse);
            reject(new Error('Socket request timeout'));
        }, time) : undefined;

        const onChatResponse = (data: any) => {
            clearTimeout(timeout);
            try {
                const parsed = completedChatResponseSchema.parse(data);
                socket.off('chat_response', onChatResponse);
                resolve(parsed);
            } catch (error) {
                reject(error);
            }
        };

        const onError = (error: any) => {
            clearTimeout(timeout);
            socket.off('chat_response', onChatResponse);
            socket.off('chat_error', onError);
            reject(error);
        };

        socket.on('chat_response', onChatResponse);
        socket.on('chat_error', onError);
    });

    return response;
}
// Parsing view name from the message's code
// todo: Get view name from chat_log.json or use some other logic
const parseViewName = (message: string) => {
    const match = /view_name\s*=\s*"([^"]+)"/.exec(message);
    if (match) {
        return match[1];
    }
    return undefined;
};

// Create conversation entry
const createConversation = (
        convMap: ConversationMap, 
        conversationId: string, 
        query: string
    ): ConversationLog => {
    if (!convMap[conversationId]) {
        const isLegacy = conversationId === 'legacy';
        convMap[conversationId] = {
            logText: isLegacy ? `(LEGACY) ${query}` : query,
            logLength: 0,
            messages: [],
        };
    }
    return convMap[conversationId];
};

// Create message pair
const createMessagePair = (log: ChatLogItem, conversationId: string) => {
    const id = generateId();
    const viewName = log?.view_name || parseViewName(log.response);
    
    const userMessage: ChatMessage = {
        conversationId,
        id,
        sender: 'user',
        text: log.query,
    };
    
    const botMessage: ChatMessage = {
        conversationId,
        id,
        sender: 'bot',
        text: log.response,
        view: viewName,
        error: log?.error,
    };
    
    return [userMessage, botMessage];
};

// Mock suggested questions
const SUGGESTED_QUESTIONS = [
    // "What columns are available in this project?",
    "Show me a summary of the cell data.",
    "Show me a scatterplot of the UMAP coordinates, colored by cluster."
    // "How many rows are in the dataset?",
    // "What is the distribution of ages?",
    // "List all unique values in the 'region' column."
];

// todo: The architecture of this hook is a bit poor, there is no single source of truth, we can have both messages and chatLog which could make it out of sync
const useChat = () => {
    // { root } is problematic here, need to revise so that we have something sensible
    // -- also test with the app not being at the root of the server
    const { projectName, mainApiRoute, projectApiRoute } = useProject(); //todo add viewName to ProjectProvider
    // const route = `${projectApiRoute}chat`;
    const routeInit = `${projectApiRoute}chat_init`;
    //! still need to consider that some things really are routes while others are just names
    const chatProgressEvent = 'chat_progress';
    const chatlogEvent = 'chat';
    // use sessionStorage to remember things like whether the chat is open or closed... localStorage for preferences like theme,
    //! but this might not be such a good idea for things like chat logs with sensitive information
    // const [messages, setMessages] = useState<ChatMessage[]>(JSON.parse(sessionStorage.getItem('chatMessages') as any) || []);
    const [messages, setMessages] = useState<ChatMessage[]>([]);
    const [isSending, setIsSending] = useState<boolean>(false);
    const [isLoadingInit, setIsLoadingInit] = useState<boolean>(false);
    const [currentRequestId, setCurrentRequestId] = useState<string>("");
    const [isInit, setIsInit] = useState<boolean>(false);
    const [requestProgress, setRequestProgress] = useState<ChatProgress | null>(null);
    const [verboseProgress, setVerboseProgress] = useState([""]);
    const cm = useChartManager();
    const viewManager = useViewManager();
    const [conversationId, setConversationId] = useState<string>(generateConversationId());
    const [conversationMap, setConversationMap] = useState<ConversationMap>({});
    const [chatLog, setChatLog] = useState<ChatLogItem[]>([]);
    const [suggestedQuestions, setSuggestedQuestions] = useState<string[]>([]);
    const socket = useMemo(() => {
        if (!cm.ipc || !cm.ipc.socket) return null;
        return cm.ipc.socket;
    }, [cm.ipc]);

    // Use React Query to manage chat logs
    const { data: chatLogData = [], isLoading: isChatLogLoading, isSuccess } = useQuery({
        queryKey: ['chatLog'],
        queryFn: async () => {
            try {
                const response = await axios.get(`${projectApiRoute}chat_log.json`);
                const parsedResponse = chatLogSchema.parse(response.data);
                return parsedResponse;
            } catch (error) {
                return [];
            }
        },
        refetchInterval: 5000, // Refetch every 5 seconds
    });

    // Update chatLog state
    useEffect(() => {
        if (chatLogData && isSuccess) {
            setChatLog(chatLogData);
        }
    }, [chatLogData, isSuccess]);

    // Initialise conversation map
    useEffect(() => {
        if (chatLog.length > 0) {
            const newConversationMap: ConversationMap = {};
        
            // Process each log entry
            chatLog.forEach((log) => {
                const conversationId = log.conversation_id || 'legacy';
                const conversation = createConversation(newConversationMap, conversationId, log.query);
                const messages = createMessagePair(log, conversationId);
                
                conversation.messages.push(...messages);
                conversation.logLength++;
            });
            
            setConversationMap(newConversationMap);
        }
    }, [chatLog]);

    // Update messages when conversation id changes
    // biome-ignore lint/correctness/useExhaustiveDependencies: We only need to change messages when conversationId changes, not everytime the conversationMap changes
    useEffect(() => {
        if (conversationMap && conversationMap[conversationId]?.messages?.length > 0) {
            const newMessages: ChatMessage[] = conversationMap[conversationId]?.messages || [];
            setMessages(newMessages);
        }
    }, [conversationId]);

    // Progress Listener
    const progressListener = useCallback((data: any) => {
        try {
            const parsed = chatProgressSchema.parse(data);
            if (parsed.id === currentRequestId) {
                setRequestProgress(parsed);
            }
        } catch (error) {
            console.error('Error parsing chat progress', error);
            setRequestProgress({
                id: currentRequestId,
                progress: 0,
                delta: 100,
                message: data,
            });
        }
    }, [currentRequestId]);

    // Verbose Progress
    const handleVerboseProgress = useCallback((msg: string) => {
        setVerboseProgress(v => [...v, msg].slice(-5));
    }, []);

    // Initiate chat function
    const chatInit = useCallback(async () => {
        if (isSending || isInit || isChatLogLoading) return;
        
        try {
            setIsSending(true);
            setIsLoadingInit(true);
            const id = generateId();
            setCurrentRequestId(id);
            //! todo - we might use socket here, in which case we need to have a different routeInit
            const response = await sendChatInitHttp('', id, routeInit, conversationId);
            if (!response) return;
            if (response?.error) throw response.message;
            
            // Only set initial message if we don't have any messages yet
            if (messages.length === 0) {
                setMessages([{
                    text: response.message,
                    sender: 'system',
                    id: generateId(),
                    conversationId,
                }]);
            }
            // todo: Update when the endpoint is ready
            // const suggestedQuestions = await axios.get("");
            if (response?.suggested_questions) setSuggestedQuestions(response?.suggested_questions);
        } catch (error: any) {
            const errorMessage =  error?.message ? 
                error.message : 
                (typeof error === "string" ? error : "ERROR: An unknown error occurred. Please try again later.")
            console.error('Error sending welcome message: ', errorMessage);
            if (messages.length === 0) {
                setMessages([{
                    text: errorMessage,
                    sender: 'system',
                    id: generateId(),
                    conversationId,
                    error: true,
                }]);
            }
        } finally {
            setIsSending(false);
            setCurrentRequestId('');
            setIsInit(true);
            setIsLoadingInit(false);
        }
    }, [isSending, isInit, routeInit, conversationId, messages.length, isChatLogLoading]);

    // Socket connection and Init chat
    useEffect(() => {
        if (!socket) return;
        if (!socket.connected) {
            console.log('Socket not connected, skipping listener registration');
            return;
        }

        // todo: Add proper error handling
        const handleSocketError = (error: any) => {
            console.error('Socket error in chat:', error);
        };
        
        socket.on('error', handleSocketError);
        socket.on(chatProgressEvent, progressListener);
        socket.on(chatlogEvent, handleVerboseProgress);

        if (!isInit && !isSending) {
            chatInit();
        }

        return () => {
            socket?.off('error', handleSocketError);
            socket?.off(chatProgressEvent, progressListener);
            socket?.off(chatlogEvent, handleVerboseProgress);
        };
    }, [socket, progressListener, handleVerboseProgress, chatInit, isInit, isSending]);
    
    // Send Message API
    const sendAPI = useCallback(async (input: string) => {
        if (!input.trim()) return;
        if (!socket) return;
        const id = generateId();
        const viewName = getViewName();
        console.log(`sending chat '${input}' from '${viewName}'`);

        try {
            setIsSending(true);
            setCurrentRequestId(id);
            setMessages(prev => [...prev, {
                text: input,
                sender: 'user',
                id: generateId(),
                conversationId,
            }])
            
            const response = await sendMessageSocket(input, id, "", conversationId);

            if (response) {
                const allViews = viewManager.all_views;
                // Add new view if it doesn't exist
                const viewName = response?.view ?? parseViewName(response.message);
                if (viewName && !allViews.includes(viewName)) {
                    viewManager.setAllViews([...allViews, viewName])
                }
                setMessages(prev => [...prev, {
                    text: response.message,
                    sender: 'bot',
                    id: generateId(),
                    conversationId,
                    view: response?.view,
                }])
            }
        } catch (error: any) {
            const errorMessage =  error?.message ? 
                error.message : 
                (typeof error === "string" ? error : "ERROR: An unknown error occurred. Please try again later.")
            setMessages(prev => [...prev, {
                text: errorMessage,
                sender: 'bot',
                id: generateId(),
                conversationId,
                error: true,
            }]);
            console.error("Error sending message: ", errorMessage);
        } finally {
            // queryClient.invalidateQueries({ queryKey: ['chatLog'] });
            setIsSending(false);
            setCurrentRequestId('');
            setRequestProgress(null);
        }
    }, [conversationId, socket, viewManager]);


    // Start New Conversation
    const startNewConversation = useCallback(async () => {
        const newConversationId = generateConversationId();
        const id = generateId();
        try {
            setConversationId(newConversationId);
            setIsLoadingInit(true);
            setCurrentRequestId(id);
            setIsSending(true);
            const response = await sendChatInitHttp('', id, routeInit, newConversationId);
            if (!response) return;

            if (response?.error) throw response.message;
            // Only set initial message if we don't have any messages yet
            setMessages([{
                text: response.message,
                sender: 'system',
                id: generateId(),
                conversationId: newConversationId
            }]);
        } catch (error: any) {
            const errorMessage =  error?.message ? 
                error.message : 
                (typeof error === "string" ? error : "ERROR: An unknown error occurred. Please try again later.")
            setMessages([{
                text: errorMessage,
                sender: 'system',
                id: generateId(),
                conversationId: newConversationId,
                error: true,
            }]);
            console.error("Error starting new conversation: ", errorMessage);
        } finally {
            setIsLoadingInit(false);
            setCurrentRequestId("");
            setIsSending(false);
        }
    }, [routeInit]);

    // Switch conversation
    const switchConversation = useCallback((id: string) => {
        console.log(`switching conversation to ${id}`);
        setCurrentRequestId("");
        setIsSending(false);
        setRequestProgress(null);
        setConversationId(id);
    }, []);

    return {
        messages,
        isSending,
        sendAPI,
        requestProgress,
        verboseProgress,
        isChatLogLoading,
        conversationId,
        startNewConversation,
        switchConversation,
        conversationMap,
        isLoadingInit,
        suggestedQuestions,
    };
};

/**
 * This function is used to navigate to a different view.
 * 
 * @param view - The name of the view to navigate to.
 * @param needsRefresh - Whether the page will need to be refreshed after navigating - for example,
 * if datasources have changed. We should be able to determine this automatically in the future
 * (and ideally not need to refresh the page).
 * 
 * ! doesn't belong in this file, I intend to refactor soon, but for quick prototyping it's here
 */
export async function navigateToView(view: string, needsRefresh = false, callback?: () => void) {
    if (!needsRefresh) {
        window.mdv.chartManager.changeView(view);
    } else {
        const params = new URLSearchParams(window.location.search);
        params.set("view", view);
        window.history.replaceState(
            {},
            "",
            `${window.location.pathname}?${params}`,
        );
        window.location.reload();
    }
    callback?.();
}

export default useChat;