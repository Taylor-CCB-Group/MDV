import { useProject } from "@/modules/ProjectContext";
import axios from "axios";
import { useCallback, useEffect, useState, useRef } from "react";
import { z } from 'zod';
import _ from 'lodash';
import { useQuery, useQueryClient } from '@tanstack/react-query';

const completedChatResponseSchema = z.object({
    message: z.string(),
    /** 
     * This should be a string that identifies newly created views.
     */
    view: z.optional(z.string()).describe('The name of the newly created view, if a view was created'),
    // timestamp: z.string(),
});

const chatProgressSchema = z.object({
    id: z.string().describe('The ID of the message this progress is associated with'),
    progress: z.number(),
    delta: z.number(),
    message: z.string(),
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
    conversation_id: z.optional(z.string()),
    timestamp: z.optional(z.string()),
});
const chatLogSchema = z.array(chatLogItemSchema);

export type ChatLogItem = z.infer<typeof chatLogItemSchema>;

export const useChatLog = () => {
    const { root } = useProject();
    const route = `${root}/chat_log.json`;

    const { data: chatLog = [], isLoading } = useQuery({
        queryKey: ['chatLog'],
        queryFn: async () => {
            const response = await axios.get(route);
            return chatLogSchema.parse(response.data);
        },
        refetchInterval: 5000,
    });

    return { chatLog, isLoading };
}

type ChatResponse = z.infer<typeof completedChatResponseSchema>;

export type ChatMessage = {
    text: string;
    view?: string;
    sender: 'user' | 'bot' | 'system';
    id: string;
    conversationId: string;
};

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


const sendMessage = async (message: string, id: string, route = '/chat', conversationId?: string) => {
    // we should send information about the context - in particular, which view we're in
    // could consider a streaming response here rather than socket
    const response = await axios.post<ChatResponse>(route, { message, id, conversation_id: conversationId });
    const parsed = completedChatResponseSchema.parse(response.data); // may throw an error if the response is not valid
    return parsed;
};

const DefaultMessage: ChatMessage = {
    text: 'Hello! How can I help you?',
    sender: 'bot',
    id: generateId(),
    conversationId: generateConversationId()
}

export type ConversationLog = {
    logText: string,
    logLength: number,
    messages: ChatMessage[],
}

export type ConversationMap = Record<string, ConversationLog>;

const useChat = () => {
    // { root } is problematic here, need to revise so that we have something sensible
    // -- also test with the app not being at the root of the server
    const { projectName, mainApiRoute, projectApiRoute } = useProject(); //todo add viewName to ProjectProvider
    const route = `${projectApiRoute}chat`;
    const routeInit = `${projectApiRoute}chat_init`;
    //! still need to consider that some things really are routes while others are just names
    const progressRoute = `/project/${projectName}/chat_progress`;
    const verboseRoute = `/project/${projectName}/chat`;
    // use sessionStorage to remember things like whether the chat is open or closed... localStorage for preferences like theme,
    //! but this might not be such a good idea for things like chat logs with sensitive information
    // const [messages, setMessages] = useState<ChatMessage[]>(JSON.parse(sessionStorage.getItem('chatMessages') as any) || []);
    const [messages, setMessages] = useState<ChatMessage[]>([]);
    const [isSending, setIsSending] = useState<boolean>(false);
    const [currentRequestId, setCurrentRequestId] = useState<string>("");
    const [isInit, setIsInit] = useState<boolean>(false);
    const [requestProgress, setRequestProgress] = useState<ChatProgress | null>(null);
    const [verboseProgress, setVerboseProgress] = useState([""]);
    const cm = window.mdv.chartManager;
    const [conversationId, setConversationId] = useState<string>(generateConversationId());
    const [conversationMap, setConversationMap] = useState<ConversationMap>({});
    const [chatLog, setChatLog] = useState<ChatLogItem[]>([]);
    const queryClient = useQueryClient();

    // Use React Query to manage chat logs
    const { data: chatLogData = [], isLoading: isChatLogLoading, isSuccess } = useQuery({
        queryKey: ['chatLog'],
        queryFn: async () => {
            const response = await axios.get(`${projectApiRoute}chat_log.json`);
            const parsedResponse = chatLogSchema.parse(response.data);
            return parsedResponse;
        },
        refetchInterval: 5000, // Refetch every 5 seconds
    });

    useEffect(() => {
        if (!chatLogData || !isSuccess) return;
        
        // Create a set of stringified items for comparison
        const currentItems = new Set(chatLog?.map(item => JSON.stringify(item)) || []);
        const newItems = new Set(chatLogData.map(item => JSON.stringify(item)));
        
        // Check if the sets have the same size and all items from newItems are in currentItems
        const hasChanges = chatLog?.length !== chatLogData.length || 
            !Array.from(newItems).every(item => currentItems.has(item));
        
        if (hasChanges) {
            setChatLog(chatLogData);
        }
    }, [chatLogData, chatLog, isSuccess]);

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

    const handleVerboseProgress = useCallback((msg: string) => {
        setVerboseProgress(v => [...v, msg].slice(-5));
    }, []);

    const chatInit = useCallback(async () => {
        if (isSending || isInit) return;
        
        setIsSending(true);
        try {
            const id = generateId();
            setCurrentRequestId(id);
            const response = await sendMessage('', id, routeInit, conversationId);
            // Only set initial message if we don't have any messages yet
            if (messages.length === 0) {
                setMessages([{
                    text: response.message,
                    sender: 'system',
                    id: generateId(),
                    conversationId
                }]);
            }
        } catch (error) {
            console.error('Error sending welcome message', error);
        }
        setIsSending(false);
        setCurrentRequestId('');
        setIsInit(true);
    }, [isSending, isInit, routeInit, conversationId, messages.length]);

    useEffect(() => {
        if (!cm.ipc) return;
        
        const { socket } = cm.ipc;
        socket?.on(progressRoute, progressListener);
        socket?.on(verboseRoute, handleVerboseProgress);

        if (!isInit && !isSending) {
            chatInit();
        }

        return () => {
            socket?.off(progressRoute, progressListener);
            socket?.off(verboseRoute, handleVerboseProgress);
        };
    }, [cm.ipc, progressRoute, verboseRoute, progressListener, handleVerboseProgress, chatInit, isInit, isSending]);

    
    useEffect(() => {
        if (chatLog.length > 0) {
            const newConversationMap: ConversationMap = {};
            const newChatLog = chatLog;
            newChatLog.forEach((log) => {
                const id = generateId();
                if (log.conversation_id) {

                    if (!newConversationMap[log.conversation_id]) {
                        newConversationMap[log.conversation_id] = {
                            logText: log.query,
                            logLength: 0,
                            messages: [],
                        };
                    }
                    const userMessage: ChatMessage = {
                        conversationId: log.conversation_id,
                        id,
                        sender: 'user',
                        text: log.query,
                    }

                    const botMessage: ChatMessage = {
                        conversationId: log.conversation_id,
                        id,
                        sender: 'bot',
                        text: `I ran some code for you:\n\n\`\`\`python\n${log.response}\n\`\`\``,
                    }

                    newConversationMap[log.conversation_id].messages.push(...[userMessage, botMessage]);
                    newConversationMap[log.conversation_id].logLength++;
                } else {
                    if (!newConversationMap['legacy']) {
                        newConversationMap['legacy'] = {
                            logText: `(LEGACY) ${log.query}`,
                            logLength: 0,
                            messages: [],
                        };
                    }
                    
                    const userMessage: ChatMessage = {
                        conversationId: 'legacy',
                        id,
                        sender: 'user',
                        text: log.query,
                    }
                    
                    const botMessage: ChatMessage = {
                        conversationId: 'legacy',
                        id,
                        sender: 'bot',
                        text: `I ran some code for you:\n\n\`\`\`python\n${log.response}\n\`\`\``,
                    }
                    
                    newConversationMap['legacy'].messages.push(...[userMessage, botMessage]);
                    newConversationMap['legacy'].logLength++;
                }
            })
            setConversationMap(newConversationMap);
        }
    }, [chatLog]);

    useEffect(() => {
        const newMessages: ChatMessage[] = conversationMap[conversationId]?.messages || [];
        setMessages(newMessages);
    }, [conversationId, conversationMap]);
    
    const sendAPI = useCallback(async (input: string) => {
        if (!input.trim()) return;
        
        const id = generateId();
        const viewName = getViewName();
        console.log(`sending chat '${input}' from '${viewName}'`);

        try {
            setIsSending(true);
            setCurrentRequestId(id);
            await sendMessage(input, id, route, conversationId);
            queryClient.invalidateQueries({ queryKey: ['chatLog'] });
        } catch (error) {
            setMessages(prev => [...prev, {
                text: `Error: ${error}`,
                sender: 'bot',
                id,
                conversationId
            }]);
        }
        setCurrentRequestId('');
        setRequestProgress(null);
        setIsSending(false);
    }, [conversationId, route, queryClient]);

    const startNewConversation = useCallback(async () => {
        const newConversationId = generateConversationId();
        setConversationId(newConversationId);
        try {
            const response = await sendMessage('', generateId(), routeInit, newConversationId);
            // Only set initial message if we don't have any messages yet
            setMessages([{
                text: response.message,
                sender: 'system',
                id: generateId(),
                conversationId: newConversationId
            }]);
        } catch (error) {
            console.log("error: ", error);
        }
    }, [routeInit]);

    const switchConversation = useCallback((id: string) => {
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
export async function navigateToView(view: string, needsRefresh = false) {
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
}

export default useChat;