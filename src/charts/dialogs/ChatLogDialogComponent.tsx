import { useState } from "react";
import { type ChatLogItem, useChatLog } from "./ChatAPI";
import ReactMarkdown from "react-markdown";
import SyntaxHighlighter from 'react-syntax-highlighter';
import { dracula } from 'react-syntax-highlighter/dist/esm/styles/hljs';
import {
    Accordion,
    AccordionContent,
    AccordionItem,
    AccordionTrigger,
} from "@/components/ui/accordion"


const Response = ({ response }: { response: string }) => {
    const [expanded, setExpanded] = useState(false);
    return (
        <div onClick={() => setExpanded(!expanded)}
            className={`${expanded ? '' : 'max-h-12 overflow-hidden'}`}
        >
            {/* {response} */}
            {/* biome-ignore lint/correctness/noChildrenProp: <explanation> */}
            <ReactMarkdown children={response} />
        </div>
    )
}

const Code = ({ code }: { code: string }) => {
    const [expanded, setExpanded] = useState(false);
    return (
        <div onClick={() => setExpanded(!expanded)}
            className={`${expanded ? '' : 'max-h-12 overflow-hidden'}`}
        >
            {/* biome-ignore lint/correctness/noChildrenProp: <explanation> */}
            <SyntaxHighlighter language={"python"} children={code} style={dracula} />
        </div>
    )
}

const LogTableItem = ({ item }: { item: ChatLogItem }) => {
    return (
        <tr>
            <td><Response response={item.context} /></td>
            <td><Response response={item.query} /></td>
            <td><Response response={item.prompt_template} /></td>
            <td><Code code={item.response} /></td>
        </tr>
    )
}

const ChatLogTable = () => {
    const { chatLog } = useChatLog();
    return (
        <div>
            <table>
                <thead>
                    <tr>
                        <th>Context</th>
                        <th>Query</th>
                        <th>Prompt template</th>
                        <th>Response</th>
                    </tr>
                </thead>
                <tbody>
                    {chatLog.map((item, index) => (
                        <LogTableItem key={index} item={item} />
                    ))}
                </tbody>
            </table>
        </div>
    )
}

const LogItem = ({ item }: { item: ChatLogItem }) => {
    return (
        <AccordionItem className="w-full" value={item.query}>
            <AccordionTrigger>{item.query}</AccordionTrigger>
            <AccordionContent>
                <Code code={item.response} />
            </AccordionContent>
        </AccordionItem>
    )
}

const ChatLog = () => {
    const { chatLog } = useChatLog();
    return (
        <Accordion type='multiple' collapsible>
            {chatLog.map((item) => (
                <LogItem key={item.query} item={item} />
            ))}
        </Accordion>
    )
}


export default ChatLog;
