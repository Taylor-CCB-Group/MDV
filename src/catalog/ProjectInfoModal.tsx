import {
    Box,
    Dialog,
    DialogContent,
    DialogTitle,
} from "@mui/material";
import type React from "react";
import { DialogCloseIconButton } from "./ProjectRenameModal";
import ReactMarkdown from "react-markdown";
import remarkGfm from "remark-gfm";
import rehypeRaw from "rehype-raw";
import rehypeSanitize from "rehype-sanitize";
import remarkBreaks from "remark-breaks";

export interface ProjectInfoModalProps {
    open: boolean;
    onClose: () => void;
    readme?: string;
}

const ProjectInfoModal: React.FC<ProjectInfoModalProps> = ({
    open,
    onClose,
    readme,
}) => {
    return (
        <Dialog
            open={open}
            onClose={onClose}
            aria-labelledby="project_info_dialog-title"
            maxWidth="sm"
            fullWidth
        >
            <DialogTitle id="project_info_dialog-title">
                Project Information
                <DialogCloseIconButton onClose={onClose} />
            </DialogTitle>
            <DialogContent dividers>
                {readme ? (
                    <Box className="readme-markdown">
                        <ReactMarkdown
                            remarkPlugins={[remarkGfm, remarkBreaks]} 
                            rehypePlugins={[rehypeRaw, rehypeSanitize]}
                            components={{
                                // Open the links in new tab
                                a: ({node, ...props}) => (
                                    <a {...props} target="_blank" rel="noopener noreferrer" />
                                )
                            }}
                        >
                            {readme}
                        </ReactMarkdown>
                    </Box>
                ) : (
                    <div>No information available.</div>
                )}
            </DialogContent>
        </Dialog>
    );
};

export default ProjectInfoModal;
