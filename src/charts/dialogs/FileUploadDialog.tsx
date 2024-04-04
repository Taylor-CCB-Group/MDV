import React, { useState, useCallback, useReducer } from "react";
import styled from "styled-components";
import { useDropzone } from "react-dropzone";
import { createRoot } from "react-dom/client";
import { BaseDialog } from "../../utilities/Dialog.js";

// Styled components
const Container = styled.div`
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  padding: 10px;
`;

const StatusContainer = styled.div`
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  width: 100%;
  height: 150px;
`;

const SuccessContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  background-color: #f0f8ff;
  box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
  border: 1px solid #e0e0e0;
`;

const SuccessHeading = styled.h1`
  color: #333;
  margin-bottom: 15px;
`;

const SuccessText = styled.p`
  font-size: 20px;
  color: #555;
  margin-bottom: 20px;
  text-align: center;
`;

const RefreshButton = styled.button`
  padding: 10px 20px;
  background-color: green;
  color: white;
  cursor: pointer;
  border-radius: 5px;
  margin-top: 20px;
  margin-bottom: 20px;
  align-self: center;
`;

const DropzoneContainer = styled.div<{ isDragOver: boolean }>`
  border: 2px dashed gray;
  padding: 10px;
  text-align: center;
  background-color: ${(props) => (props.isDragOver ? "lightgray" : "white")};
  min-width: 90%;
  min-height: 90%;
  max-height: 90%;
`;

const FileInputLabel = styled.label`
  padding: 10px 20px;
  border: 1px solid #ccc;
  border-radius: 4px;
  cursor: pointer;
  background-color: #f0f0f0;
  display: inline-block;
  margin: 10px 0;
`;

const Button = styled.button<{ color?: string; disabled?: boolean; size?: string; marginTop?: string; }>`
  align-self: center;
  padding: ${(props) => props.size || "10px 20px"};
  background-color: ${(props) => props.color || "blue"};
  color: white;
  cursor: pointer;
  margin-top: ${(props) => props.marginTop || "10px"};
  opacity: ${(props) => (props.disabled ? 0.5 : 1)};
`;

const ProgressBar = styled.progress`
  width: 100%;
  height: 30px;
  margin-bottom: 20px;
  background-color: #f3f3f3;
  border: 1px solid #ccc;
  border-radius: 5px;
`;

const Message = styled.p`
  font-size: 18px;
  font-weight: bold;
  color: #333;
  font-family: Arial, sans-serif;
  text-align: center;
`;

const FileSummary = styled.div`
  border-radius: 8px;
  background-color: #f0f8ff;
  box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
  border: 1px solid #e0e0e0;
  max-width: 400px;
  width: 90%;
  text-align: center;
`;

const FileSummaryHeading = styled.h2`
  color: #333;
  margin: 0 0 20px 0;
`;

const FileSummaryText = styled.p`
  font-size: 16px;
  margin: 10px 0;
  padding: 0;
  color: #555;
`;

const ErrorContainer = styled.div`
  background-color: #fff0f0; // Light red background for error visibility
  color: #d8000c; // Dark red text color for contrast and readability
  border: 1px solid #ffbaba; // Lighter red border color
  padding: 20px;
  margin: 20px 0;
  border-radius: 5px;
  box-shadow: 0 2px 4px rgba(0,0,0,0.1); // Slight shadow for depth
  font-family: Arial, sans-serif;
  font-size: 14px;
  text-align: left;
  width: 90%; // Responsive width
  max-width: 600px; // Max width to avoid too wide error messages
`;

const ErrorHeading = styled.h4`
  margin-top: 0;
  font-size: 16px;
  font-weight: bold;
`;

// Reducer function
const reducer = (state, action) => {
  switch (action.type) {
    case "SET_SELECTED_FILES":
      return { ...state, selectedFiles: action.payload };
    case "SET_IS_UPLOADING":
      return { ...state, isUploading: action.payload };
    case "SET_FILE_SUMMARY":
      return { ...state, fileSummary: action.payload };
    case "SET_IS_INSERTING":
      return { ...state, isInserting: action.payload };
    case "SET_SUCCESS":
      return { ...state, success: action.payload };
    case "SET_ERROR":
      return { ...state, error: action.payload };
    default:
      return state;
  }
};

// Constants
const UPLOAD_INTERVAL = 30;
const UPLOAD_DURATION = 3000;
const UPLOAD_STEP = 100 / (UPLOAD_DURATION / UPLOAD_INTERVAL);

interface FileUploadDialogComponentProps {
  onClose: () => void;
}

// Custom hook for file upload progress
const useFileUploadProgress = () => {
  const [progress, setProgress] = useState(0);

  const startProgress = useCallback(() => {
    setProgress(0);
    const interval = setInterval(() => {
      setProgress((prevProgress) => {
        const updatedProgress = prevProgress + UPLOAD_STEP;
        if (updatedProgress >= 100) {
          clearInterval(interval);
          return 100;
        }
        return updatedProgress;
      });
    }, UPLOAD_INTERVAL);
  }, []);

  const resetProgress = useCallback(() => {
    setProgress(0);
  }, []);

  return { progress, startProgress, resetProgress };
};

const FileUploadDialogComponent: React.FC<FileUploadDialogComponentProps> = ({ onClose }) => {
  const [state, dispatch] = useReducer(reducer, {
    selectedFiles: [],
    isUploading: false,
    fileSummary: null,
    isInserting: false,
    success: false,
    error: null,
  });

  const { progress, startProgress, resetProgress } = useFileUploadProgress();

  const onDrop = useCallback((acceptedFiles: File[]) => {
    dispatch({ type: "SET_SELECTED_FILES", payload: acceptedFiles });
  }, []);

  const { getRootProps, getInputProps, isDragActive } = useDropzone({ onDrop });

  const handleUploadClick = async () => {
    if (!state.selectedFiles.length) {
      dispatch({ type: "SET_ERROR", payload: "noFilesSelected" });
      return;
    }

    dispatch({ type: "SET_IS_UPLOADING", payload: true });
    resetProgress();
    dispatch({
      type: "SET_ERROR",
      payload: {
        message: "Upload failed due to a network error.",
        traceback: "Error: Network Error at uploadFileFunction"
      }
    });

    startProgress();

    try {
      // Simulating file upload

      // const endpoint = 'https://example.com/upload';
      // const formData = new FormData();

      // Array.from(selectedFiles).forEach(file => {
      //   formData.append('files', file);
      // });

      // try {
      //   const response = await fetch(endpoint, {
      //     method: 'POST',
      //     body: formData,
      //   });

      //   if (response.ok) {
      //     const result = await response.json();
      //     console.log("Upload successful:", result);
      //   } else {
      //     console.error("Upload failed:", response.statusText);
      //   }
      // } catch (error) {
      //   console.error("Error during upload:", error);
      // } finally {
      //   setIsUploading(false); // Finish uploading
      // }

      await new Promise((resolve) => setTimeout(resolve, UPLOAD_DURATION));

      dispatch({ type: "SET_IS_UPLOADING", payload: false });
      dispatch({
        type: "SET_FILE_SUMMARY",
        payload: {
          rows: 124,
          columns: 12,
          size: (10241024 / (1024 * 1024)).toFixed(2),
        },
      });
    } catch (error) {
      dispatch({ type: "SET_IS_UPLOADING", payload: false });
      dispatch({
        type: "SET_ERROR",
        payload: {
          message: "Upload failed due to a network error.",
          traceback: "Error: Network Error at uploadFileFunction"
          // The traceback here is just an example; it could be the actual error stack or details you get from your upload logic
        }
      });
    }
  };

  const handleInsertClick = async () => {
    dispatch({ type: "SET_FILE_SUMMARY", payload: null });
    dispatch({ type: "SET_IS_INSERTING", payload: true });
    resetProgress();
    dispatch({
      type: "SET_ERROR",
      payload: {
        message: "Upload failed due to a network error.",
        traceback: "Error: Network Error at uploadFileFunction"
        // The traceback here is just an example; it could be the actual error stack or details you get from your upload logic
      }
    });

    startProgress();

    // const endpoint = 'https://example.com/upload';
    // const formData = new FormData();

    // Array.from(selectedFiles).forEach(file => {
    //   formData.append('files', file);
    // });

    // try {
    //   const response = await fetch(endpoint, {
    //     method: 'POST',
    //     body: formData,
    //   });

    //   if (response.ok) {
    //     const result = await response.json();
    //     console.log("Upload successful:", result);
    //   } else {
    //     console.error("Upload failed:", response.statusText);
    //   }
    // } catch (error) {
    //   console.error("Error during upload:", error);
    // } finally {
    //   setIsUploading(false); // Finish uploading
    // }

    try {
      // Simulating file insertion
      await new Promise((resolve) => setTimeout(resolve, UPLOAD_DURATION));

      dispatch({ type: "SET_IS_INSERTING", payload: false });
      dispatch({ type: "SET_SUCCESS", payload: true });
    } catch (error) {
      dispatch({ type: "SET_IS_INSERTING", payload: false });
      dispatch({ type: "SET_ERROR", payload: "insertError" });
    }
  };

  return (
    <Container>
      {state.isUploading ? (
        <StatusContainer>
          <Message>{"Your file is being uploaded, please wait..."}</Message>
          <ProgressBar value={progress} max="100" />
        </StatusContainer>
      ) : state.fileSummary ? (
        <>
          <FileSummary>
            <FileSummaryHeading>{"Uploaded File Summary"}</FileSummaryHeading>
            <FileSummaryText>
              <strong>{"File name"}</strong> {state.selectedFiles[0].name}
            </FileSummaryText>
            <FileSummaryText>
              <strong>{"Number of rows"}</strong> {state.fileSummary.rows}
            </FileSummaryText>
            <FileSummaryText>
              <strong>{"Number of columns"}</strong> {state.fileSummary.columns}
            </FileSummaryText>
            <FileSummaryText>
              <strong>{"File size"}</strong> {state.fileSummary.size} MB
            </FileSummaryText>
          </FileSummary>
          <div
            style={{
              display: 'flex',
              justifyContent: 'center',
              alignItems: 'center',
              gap: '60px',
            }}
          >
            <Button onClick={handleInsertClick}>{"Confirm"}</Button>
            <Button color="red" size="10px 30px" onClick={onClose}>
              {"Cancel"}
            </Button>
          </div>
        </>
      ) : state.isInserting ? (
        <div
          style={{
            display: 'flex',
            flexDirection: 'column', // Stack items vertically
            justifyContent: 'center',
            alignItems: 'center',
            width: '100%', // Ensure the wrapper spans the entire width of its parent
            height: '150px', // Adjust this height as needed to accommodate both the progress bar and the message
          }}
        >
          <Message>{"Your file is being processed, please wait..."}</Message>
          <ProgressBar value={progress} max="100" />
        </div>
      ) : state.success ? (
        <>
          <SuccessContainer>
            <SuccessHeading>Success!</SuccessHeading>
            <SuccessText>
              The file was uploaded successfully to the database.
            </SuccessText>
            <RefreshButton onClick={() => window.location.reload()}>
              Refresh Page
            </RefreshButton>
          </SuccessContainer>
        </>

      ) : state.error ? (
        <>
          <ErrorContainer>
            <ErrorHeading>An error occurred while uploading the file:</ErrorHeading>
            <p>{state.error.message}</p>
            {state.error.traceback && (
              <pre>{state.error.traceback}</pre> // Use <pre> for preformatted text like tracebacks
            )}
          </ErrorContainer>
        </>
      ) : (
        <>
          <DropzoneContainer {...getRootProps()} isDragOver={isDragActive} aria-label={"dropzoneLabel"}>
            <input {...getInputProps()} />
            {isDragActive ? (
              <h2>{"dropFilesHere"}</h2>
            ) : (
              <h2>
                {state.selectedFiles.length > 0 ? "Selected file: " + state.selectedFiles[0].name : "Drag and drop files here or click the button below to upload"}
              </h2>
            )}
            <FileInputLabel htmlFor="fileInput">{"Choose File"}</FileInputLabel>
          </DropzoneContainer>
          <Button onClick={handleUploadClick} disabled={!state.selectedFiles.length || state.isUploading} marginTop="20px">
            {"Upload File"}
          </Button>
        </>
      )}
    </Container>
  );
};

class FileUploadDialogReact extends BaseDialog {
  root: ReturnType<typeof createRoot>;

  constructor() {
    super(
      {
        title: "File Upload",
        width: 450,
        height: 295,
      },
      null
    );
    this.outer.classList.add("fileUploadDialog");
    if (this.dialog) {
      this.root = createRoot(this.dialog);
      this.root.render(<FileUploadDialogComponent onClose={() => this.close()} />);
    } else {
      console.error("Dialog element not found");
    }
  }

  close(): void {
    super.close();
    if (this.root) {
      this.root.unmount();
    }
  }
}

BaseDialog.experiment["FileUploadDialogReact"] = FileUploadDialogReact;

export default FileUploadDialogReact;