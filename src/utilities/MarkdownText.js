import {marked} from "marked";
import {sanitize} from 'dompurify';

function renderText(text){
    return sanitize(marked.parse(text));
}

export default renderText