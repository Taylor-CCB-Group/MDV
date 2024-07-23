import { createRoot } from 'react-dom/client';
import App from "./catalog";

const container = document.createElement('div');
document.body.appendChild(container);
const root = createRoot(container);
root.render(<App />);