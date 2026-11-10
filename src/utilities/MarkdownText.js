import purify from "dompurify";
import { marked } from "marked";

function renderText(text) {
    if (!text) return "";
    return purify.sanitize(marked.parse(text));
}

export default renderText;
