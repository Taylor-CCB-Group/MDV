import { useProject } from "@/modules/ProjectContext";
import axios from "axios";
import { useCallback, useEffect, useState } from "react";
import { z } from 'zod';

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
        // refresh periodically - there may be a better way (quite possible involving useQuery,
        // or socketio...).
        const interval = setInterval(fetchChatLog, 5000);
        return () => clearInterval(interval);
    }, [route]);

    return { chatLog, isLoading };
}

type ChatResponse = z.infer<typeof completedChatResponseSchema>;

export type ChatMessage = {
    text: string;
    view?: string; //maybe this type should be more assocated with ChatResponse
    sender: 'user' | 'bot' | 'system';
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


const sendMessage = async (message: string, id: string, route = '/chat') => {
    // we should send information about the context - in particular, which view we're in
    // could consider a streaming response here rather than socket
    const response = await axios.post<ChatResponse>(route, { message, id });
    const parsed = completedChatResponseSchema.parse(response.data); // may throw an error if the response is not valid
    return parsed;
};

const DefaultMessage: ChatMessage = {
    text: 'Hello! How can I help you?',
    sender: 'bot',
    id: generateId(),
}

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

    const progressListener = useCallback((data: any) => {
        console.log('chat message', data);
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

    useEffect(() => {
        const { socket } = window.mdv.chartManager.ipc;

        // event-name like 'chat', room for project... id associated with original request.
        socket?.on(progressRoute, progressListener);
        const verboseProgress = (msg: string) => setVerboseProgress(v => [...v, msg].slice(-5));
        socket?.on(verboseRoute, verboseProgress);
        console.log(`addded listeners for '${progressRoute}' and '${verboseRoute}'`);
        const chatInit = async () => {
            setIsSending(true);
            try {
                const id = generateId();
                setCurrentRequestId(id);
                const response = await sendMessage('', id, routeInit);
                setMessages(messages => messages.length ? messages : [{ text: response.message, sender: 'system', id: generateId() }]);
            } catch (error) {
                console.error('Error sending welcome message', error);
            }
            setIsSending(false);
            setCurrentRequestId(''); //todo review react query etc
            setIsInit(true);
        };
        if (!isSending && !isInit) chatInit();
        return () => {
            socket?.off(progressRoute, progressListener);
            socket?.off(verboseRoute, verboseProgress);
        }
    }, [isSending, isInit, routeInit, progressRoute, verboseRoute, progressListener]);
    useEffect(() => {
        sessionStorage.setItem('chatMessages', JSON.stringify(messages));
    }, [messages]);

    const appendMessage = (message: string, sender: 'bot' | 'user', view?: string) => {
        //we should be using an id passed as part of the message, not generating one here.
        //also - id as react key if we have an id shared between query and response may be a conflict
        const msg = { text: message, sender, id: generateId(), view };
        setMessages((prevMessages) => [...prevMessages, msg]);
    };

    const sendAPI = async (input: string) => {
        if (!input.trim()) return;
        const id = generateId();
        const viewName = getViewName();
        console.log(`sending chat '${input}' from '${viewName}'`)
        appendMessage(input, 'user');

        try {
            setIsSending(true);
            setCurrentRequestId(id);
            const response = await sendMessage(input, id, route);
            // we should be appending more stuff, and then rendering appropriately, 
            // with things like a button to navigate to the view if view property is present.
            appendMessage(response.message, 'bot', response.view);
            //todo - navigating via button rather than automatically.

            // if (response.view) navigateToView(response.view);
        } catch (error) {
            appendMessage(`Error: ${error}`, 'bot');
        }
        setCurrentRequestId(id);
        setRequestProgress(null);
        setIsSending(false);
    };
    return { messages, isSending, sendAPI, requestProgress, verboseProgress };
}

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