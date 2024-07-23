import {marked} from "marked";
import purify from 'dompurify';

function renderText(text){
    if (!text) return "";
    return purify.sanitize(marked.parse(text));
}

export default renderText