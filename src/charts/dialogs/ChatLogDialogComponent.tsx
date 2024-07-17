import { useState } from "react";
import { useChatLog } from "./ChatAPI";
import ReactMarkdown from "react-markdown";
import SyntaxHighlighter from 'react-syntax-highlighter';
import { dracula } from 'react-syntax-highlighter/dist/esm/styles/hljs';


const Response = ({ response }: { response: string }) => {
    const [expanded, setExpanded] = useState(false);
    return (
        <div onClick={() => setExpanded(!expanded)}
        className={`${expanded ? '' : 'max-h-12 overflow-hidden'}`}
        >
            {/* {response} */}
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
        <SyntaxHighlighter language={"python"} children={code} style={dracula}/>
        </div>
    )
}


const ChatLog = () => {
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
                        <tr key={index}>
                            <td><Response response={item.context} /></td>
                            <td><Response response={item.query} /></td>
                            <td><Response response={item.prompt_template} /></td>
                            <td><Code code={item.response} /></td>
                        </tr>
                    ))}
                </tbody>
            </table>
        </div>
    )
}

export default ChatLog;