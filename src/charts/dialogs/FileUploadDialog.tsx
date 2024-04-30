import React, { useState, useCallback, useReducer, type PropsWithChildren } from "react";
import { useDropzone } from "react-dropzone";
import axios from 'axios';
import { useProject } from "../../modules/ProjectContext";

const Container = ({ children }: PropsWithChildren) => {
  return (
    <div className="flex flex-col content-center items-center h-max dark:bg-black-800 dark:text-white">
      {children}
    </div>
  );
};

const StatusContainer = ({ children }: PropsWithChildren) => {
  return (
    <div className="flex flex-col justify-center items-center w-full h-150 dark:bg-#333">
      {children}
    </div>
  );
};

const SuccessContainer = ({ children }) => (
  <div className="flex flex-col items-center justify-center bg-[#f0f8ff] shadow-md border border-[#e0e0e0] m-4 dark:bg-black dark:border-gray-600">
    {children}
  </div>
);

const SuccessHeading = ({ children }) => (
  <h1 className="text-[#333] mb-1 dark:text-white">{children}</h1>
);

const SuccessText = ({ children }) => (
  <p className="text-2xl text-[#555] mb-3 text-center dark:text-gray-300">{children}</p>
);

const DropzoneContainer = ({ isDragOver, children, ...props }) => (
  <div
    {...props}
    style={{ height: "90%" }}
    className={`p-2.5 z-50 text-center border-2 border-dashed ${isDragOver ? 'bg-gray-300 dark:bg-slate-800' : 'bg-white dark:bg-black'
      } min-w-[90%]`}
  >
    {children}
  </div>
);

const FileInputLabel = ({ children, ...props }) => (
  <label {...props} className="mt-8 px-5 py-2.5 border bg-stone-200 hover:bg-stone-300 rounded cursor-pointer inline-block my-2.5 dark:bg-stone-600 dark:hover:bg-stone-500">
    {children}
  </label>
);

const Spinner = () => {
  return <div className="w-16 h-16 border-8 mt-10 border-blue-500 border-dashed rounded-full animate-spin" style={{
    borderColor: "blue transparent blue transparent"
  }}></div>;
};

const colorStyles = {
  blue: {
    bgColor: "bg-blue-600",
    hoverColor: "hover:bg-blue-700",
    darkBgColor: "dark:bg-blue-800",
    darkHoverColor: "dark:bg-blue-900",
  },
  red: {
    bgColor: "bg-red-600",
    hoverColor: "hover:bg-red-700",
    darkBgColor: "dark:bg-red-800",
    darkHoverColor: "dark:bg-red-900",
  },
  green: {
    bgColor: "bg-green-600",
    hoverColor: "hover:bg-green-700",
    darkBgColor: "dark:bg-green-800",
    darkHoverColor: "dark:bg-green-900",
  },
};

const Button = ({
  onClick,
  color = "blue",
  disabled = false,
  size = "px-5 py-2.5",
  marginTop = "mt-2.5",
  children
}) => {
  const { bgColor, hoverColor, darkBgColor, darkHoverColor } = colorStyles[color] || colorStyles.blue;

  return (
    <button
      onClick={onClick}
      className={`${size} ${marginTop} ${bgColor} ${hoverColor} ${darkBgColor} ${darkHoverColor} text-white rounded self-center cursor-pointer disabled:cursor-not-allowed disabled:bg-gray-500 disabled:opacity-70`}
      disabled={disabled}
    >
      {children}
    </button>
  );
};

const ProgressBar = ({ value, max }) => (
  <progress className="w-full h-8 mb-5 mt-10 bg-gray-200 dark:bg-white-200 border border-gray-300 rounded" value={value} max={max}></progress>
);

const Message = ({ children }) => (
  <p className="text-lg font-bold text-[#333] dark:text-white text-center">{children}</p>
);

const FileSummary = ({ children }) => (
  <div className="rounded-lg bg-[#f0f8ff] dark:bg-stone-800 shadow-md border border-[#e0e0e0] max-w-[400px] w-[90%] text-center">
    {children}
  </div>
);

const FileSummaryHeading = ({ children }) => (
  <h2 className="text-gray-800 dark:text-white mb-3">{children}</h2>
)

const FileSummaryText = ({ children }) => (
  <p className="text-base text-gray-700 dark:text-white my-1">{children}</p>
);

const ErrorContainer = ({ children }) => (
  <div className="bg-[#fff0f0] text-[#d8000c] border border-[#ffbaba] p-5 my-5 rounded shadow-sm text-left w-[90%] max-w-[600px]">
    {children}
  </div>
);

const DynamicText = ({ text, className = "" }) => (
  <div className="w-96 h-20 overflow-hidden flex items-center justify-center">
    <p className={`text-center m-0 font-bold text-sm sm:text-lg md:text-xl ${className}`}>
      {text}
    </p>
  </div>
);

const ErrorHeading = ({ children }) => (
  <h4 className="mt-0 text-lg font-bold">
    {children}
  </h4>
);

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

  return { progress, setProgress, startProgress, resetProgress };
};

const FileUploadDialogComponent: React.FC<FileUploadDialogComponentProps> = ({ onClose }) => {
  const { projectName, flaskURL, isPopout } = useProject();

  const [state, dispatch] = useReducer(reducer, {
    selectedFiles: [],
    isUploading: false,
    fileSummary: null,
    isInserting: false,
    success: false,
    error: null,
  });

  const { progress, resetProgress, setProgress } = useFileUploadProgress();

  const onDrop = useCallback((acceptedFiles: File[]) => {
    dispatch({ type: "SET_SELECTED_FILES", payload: acceptedFiles });
  }, []);

  const { getRootProps, getInputProps, isDragActive, fileRejections } = useDropzone({
    onDrop,
    accept: {
      'text/csv': ['.csv']
    }
  });
  const rejectionMessage = fileRejections.length > 0 ? "Only CSV files can be selected" : "Drag and drop files here or click the button below to upload";
  const rejectionMessageStyle = fileRejections.length > 0 ? "text-red-500" : "";


  const handleUploadClick = async () => {
    console.log("Uploading file...");
    if (!state.selectedFiles.length) {
      dispatch({ type: "SET_ERROR", payload: "noFilesSelected" });
      return;
    }

    dispatch({ type: "SET_IS_UPLOADING", payload: true });
    resetProgress();

    const formData = new FormData();
    formData.append("file", state.selectedFiles[0]);
    formData.append("name", state.selectedFiles[0].name);
    formData.append("view", "TRUE");
    formData.append("backend", "TRUE");
    formData.append("replace", "TRUE");
    
    
    const config = {
      headers: {
        'Content-Type': 'multipart/form-data'
      },
      onUploadProgress: function(progressEvent) {
        const percentComplete = Math.round((progressEvent.loaded * 100) / progressEvent.total);
        setProgress(percentComplete); // Update the progress as the file uploads
      }
    };

    try {
      const response = await axios.post(`${flaskURL}${projectName}/add_datasource`, formData, config);
      console.log('File uploaded successfully', response.data);
      dispatch({ type: "SET_IS_UPLOADING", payload: false });
      
      // Extracting rows, columns, and size from response.data
      const { rows, columns, size } = response.data;
      dispatch({
        type: "SET_FILE_SUMMARY",
        payload: {
          rows: rows,
          columns: columns,
          size: size,
        },
      });
    } catch (error) {
      console.error('Error uploading file:', error);
      dispatch({ type: "SET_IS_UPLOADING", payload: false });
      dispatch({
        type: "SET_ERROR",
        payload: {
          message: "Upload failed due to a network error.",
          traceback: error.message
        }
      });
    }
  };

  const handleInsertClick = async () => {
    dispatch({ type: "SET_FILE_SUMMARY", payload: null });
    dispatch({ type: "SET_IS_INSERTING", payload: true });
  
    try {
      const response = await axios.post(`${flaskURL}${projectName}/confirm_datasource`, { confirmed: true });
      if (response.status === 200) {
        console.log('Confirmation successful', response.data);
        dispatch({ type: "SET_IS_INSERTING", payload: false });
        dispatch({ type: "SET_SUCCESS", payload: true });
      } else {
        throw new Error('Failed to confirm');
      }
    } catch (error) {
      console.error('Error during confirmation:', error);
      dispatch({ type: "SET_IS_INSERTING", payload: false });
      dispatch({
        type: "SET_ERROR", payload: {
          message: "Confirmation failed due to a network or server error.",
          traceback: error.message || "No additional error info"
        }
      });
    }
  };

  const handleClose = async () => {
    dispatch({ type: "SET_FILE_SUMMARY", payload: null });
    dispatch({ type: "SET_IS_INSERTING", payload: true });
  
    try {
      const response = await axios.post(`${flaskURL}${projectName}/confirm_datasource`, { confirmed: false });
      if (response.status === 200) {
        console.log('Cancellation successful');
        onClose(); // Call onClose to handle any local cleanup or UI changes
      } else {
        throw new Error('Failed to cancel');
      }
    } catch (error) {
      console.error('Error during cancellation:', error);
    } finally {
      dispatch({ type: "SET_IS_INSERTING", payload: false });
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
            <Button marginTop="mt-3" onClick={handleInsertClick}>{"Confirm"}</Button>
            <Button color="red" size="px-6 py-2.5" marginTop="mt-3" onClick={handleClose}>{"Cancel"}</Button>
          </div>
        </>
      ) : state.isInserting ? (
        <StatusContainer>
          <Message>{"Your file is being processed, please wait..."}</Message>
          <Spinner />
        </StatusContainer>
      ) : state.success ? (
        <>
          <SuccessContainer>
            <SuccessHeading>Success!</SuccessHeading>
            <SuccessText>
              The file was uploaded successfully to the database.
            </SuccessText>
          </SuccessContainer>
          <Button color="green" onClick={() => window.location.reload()}>
            Refresh Page
          </Button>
        </>

      ) : state.error ? (
        <>
          <ErrorContainer>
            <ErrorHeading>An error occurred while uploading the file:</ErrorHeading>
            <p>{state.error.message}</p>
            {state.error.traceback && (
              <pre>{state.error.traceback}</pre>
            )}
          </ErrorContainer>
        </>
      ) : (
        <>
          <DropzoneContainer {...getRootProps()} isDragOver={isDragActive} aria-label={"dropzoneLabel"}>
            <input {...getInputProps()} />
            {isDragActive ? (
              <DynamicText text={"Drop files here..."}/>
            ) : (
              <DynamicText text={state.selectedFiles.length > 0 ? "Selected file: " + state.selectedFiles[0].name : rejectionMessage} className={rejectionMessageStyle} />
            )}
            <FileInputLabel htmlFor="fileInput">{"Choose File"}</FileInputLabel>
          </DropzoneContainer>
          <Button onClick={handleUploadClick}
            disabled={(state.selectedFiles.length === 0 || state.isUploading)}
            marginTop="mt-5"
            color="blue"
          >
            {"Upload File"}
          </Button>
        </>
      )}
    </Container>
  );
};

export default FileUploadDialogComponent;