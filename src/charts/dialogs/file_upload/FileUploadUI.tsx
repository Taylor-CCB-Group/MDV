import {
    type PropsWithChildren,
    forwardRef,
} from "react";

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
    gray: {
        bgColor: "bg-gray-600",
        hoverColor: "hover:bg-gray-700",
        darkBgColor: "dark:bg-gray-800",
        darkHoverColor: "dark:bg-gray-900",
    },
};

type ButtonProps = {
    onClick: () => void;
    color?: "blue" | "red" | "green" | "gray";
    disabled?: boolean;
    size?: string;
    marginTop?: string;
} & PropsWithChildren;

export const Container = ({ children }: PropsWithChildren) => {
    return (
        <div className="flex flex-col content-center items-center h-max dark:text-white">
            {children}
        </div>
    );
};

export const StatusContainer = ({ children }: PropsWithChildren) => {
    return (
        <div className="flex flex-col justify-center items-center w-full h-[150px]">
            {children}
        </div>
    );
};

export const SuccessContainer = ({ children }: PropsWithChildren) => (
    <div className="flex flex-col items-center justify-center bg-[#f0f8ff] shadow-md border border-[#e0e0e0] m-4 dark:bg-black dark:border-gray-600">
        {children}
    </div>
);

export const SuccessHeading = ({ children }: PropsWithChildren) => (
    <h1 className="text-[#333] mb-1 dark:text-white">{children}</h1>
);

export const SuccessText = ({ children }: PropsWithChildren) => (
    <p className="text-2xl text-[#555] mb-3 text-center dark:text-gray-300">
        {children}
    </p>
);

export const DropzoneContainer = forwardRef(
    ({ isDragOver, children, ...props }: any, ref) => (
        <div
            {...props}
            ref={ref}
            className={`p-4 mt-2 z-50 text-center border-2 border-dashed rounded-lg ${isDragOver ? "bg-gray-300 dark:bg-slate-800" : "bg-white dark:bg-black"} min-w-[90%]`}
        >
            {children}
        </div>
    ),
);

//! warning: className isn't being merged etc, I think it just gets overriden.
// not sure if there's a better type for `htmlFor`.
export const FileInputLabel = ({
    children,
    htmlFor,
    ...props
}: PropsWithChildren & { className: string; htmlFor: string }) => (
    <label
        htmlFor={htmlFor}
        {...props}
        className="mt-2 px-5 py-2.5 border bg-stone-200 hover:bg-stone-300 rounded cursor-pointer inline-block my-2.5 dark:bg-stone-600 dark:hover:bg-stone-500"
    >
        {children}
    </label>
);

export const Spinner = () => {
    return (
        <div
            className="w-16 h-16 border-8 mt-10 border-blue-500 border-dashed rounded-full animate-spin"
            style={{
                borderColor: "blue transparent blue transparent",
            }}
        />
    );
};

export const Button = ({
    onClick,
    color = "blue",
    disabled = false,
    size = "px-5 py-2.5",
    marginTop = "mt-2.5",
    children,
}: ButtonProps) => {
    const { bgColor, hoverColor, darkBgColor, darkHoverColor } =
        colorStyles[color] || colorStyles.blue;

    return (
        <button
            type="button"
            onClick={onClick}
            className={`${size} ${marginTop} ${bgColor} ${hoverColor} ${darkBgColor} ${darkHoverColor} text-white rounded self-center cursor-pointer disabled:cursor-not-allowed disabled:bg-gray-500 disabled:opacity-70`}
            disabled={disabled}
        >
            {children}
        </button>
    );
};

//@ts-ignore CBA
export const ProgressBar = ({ value, max }) => (
    <progress
        className="w-full h-8 mb-5 mt-10 bg-gray-200 dark:bg-white-200 border border-gray-300 rounded"
        value={value}
        max={max}
    />
);

export const Message = ({ children }: PropsWithChildren) => (
    <p className="text-lg font-bold text-[#333] dark:text-white text-center">
        {children}
    </p>
);

export const FileSummary = ({ children }: PropsWithChildren) => (
    <div className="max-w-[800px] w-[100%] text-center mb-5">{children}</div>
);

export const FileSummaryHeading = ({ children }: PropsWithChildren) => (
    <h1 className="text-gray-800 dark:text-white mb-3">{children}</h1>
);

export const FileSummaryText = ({ children }: PropsWithChildren) => (
    <>
        {typeof children === "string" ? (
            <p className="text-lg text-gray-700 dark:text-white my-1">
                {children}
            </p>
        ) : (
            <div className="text-lg text-gray-700 dark:text-white my-1">
                {children}
            </div>
        )}
    </>
);

export const ErrorContainer = ({ children }: PropsWithChildren) => (
    <div className="bg-[#fff0f0] text-[#d8000c] border border-[#ffbaba] p-5 my-5 rounded shadow-sm text-left w-[90%] max-w-[600px]">
        {children}
    </div>
);

export const DynamicText = ({
    text,
    className = "",
}: { text: string; className: string }) => (
    <div className="w-96 h-20 overflow-hidden flex items-center justify-center">
        <p
            //consider using `cn()` from lib/utils
            className={`text-center m-0 font-bold text-sm sm:text-lg md:text-m ${className}`}
        >
            {text}
        </p>
    </div>
);

export const ErrorHeading = ({ children }: PropsWithChildren) => (
    <h4 className="mt-0 text-lg font-bold">{children}</h4>
);

//@ts-ignore CBA
export const DatasourceNameInput = ({ value, onChange, isDisabled }) => (
    <div className="flex-left items-center space-x-2 pr-4">
        <label className="text-lg text-gray-700 dark:text-white my-1" htmlFor="datasourceName">
            <strong>Datasouce Name:</strong>
        </label>
        <input
            id="datasourceName"
            type="text"
            value={value}
            onChange={onChange}
            className="p-2 border rounded focus:outline-none focus:ring focus:border-blue-300 dark:text-gray-300 dark:bg-gray-800 flex-grow"
            placeholder="Enter datasource name"
            disabled={isDisabled}
        />
    </div>
);
