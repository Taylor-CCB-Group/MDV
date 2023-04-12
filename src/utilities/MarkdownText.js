import {marked} from "marked";
import purify from 'dompurify';

function renderText(text){
    return purify.sanitize(marked.parse(text));
}

export default renderText