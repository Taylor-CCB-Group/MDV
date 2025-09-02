import type React from "react";
import {
    useState,
    useCallback,
    useReducer,
    type PropsWithChildren,
    forwardRef,
    useEffect,
} from "react";
import { observer } from "mobx-react-lite";
import { useProject } from "../../../modules/ProjectContext";
import {
    CloudDownload as CloudDownloadIcon,
    Folder as FolderIcon,
    Image as ImageIcon,
    TableChart as TableIcon,
    Info as InfoIcon,
} from "@mui/icons-material";
import ReusableDialog from "../ReusableDialog";
import { debugExploreWithHTTP, debugExploreZarrDataset, extractZarrMetadata } from "./zarrMetadataUtils";

// Import zarr-js - you'll need to install this: npm install zarr
// import { openGroup } from "zarr";

export const Container = ({ children }: PropsWithChildren) => {
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

const MetadataContainer = ({ children }: PropsWithChildren) => (
    <div className="flex flex-col bg-[#f8f9fa] shadow-md border border-[#e0e0e0] m-4 p-4 rounded-lg dark:bg-black dark:border-gray-600 w-full max-w-6xl">
        {children}
    </div>
);

const SectionHeader = ({ children, icon }: PropsWithChildren & { icon?: React.ReactNode }) => (
    <div className="flex items-center gap-2 mb-3 pb-2 border-b border-gray-200 dark:border-gray-600">
        {icon}
        <h3 className="text-lg font-semibold text-[#333] dark:text-white">{children}</h3>
    </div>
);

const MetadataGrid = ({ children }: PropsWithChildren) => (
    <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 mb-6">
        {children}
    </div>
);

const MetadataCard = ({ title, children }: PropsWithChildren & { title: string }) => (
    <div className="bg-white dark:bg-gray-800 border border-gray-200 dark:border-gray-600 rounded-lg p-3">
        <h4 className="font-medium text-sm text-gray-600 dark:text-gray-300 mb-2">{title}</h4>
        <div className="text-sm text-gray-800 dark:text-gray-200">{children}</div>
    </div>
);

const DataTable = ({ data, maxRows = 10 }: { data: Record<string, any>; maxRows?: number }) => {
    const entries = Object.entries(data).slice(0, maxRows);
    
    return (
        <div className="overflow-x-auto">
            <table className="min-w-full text-sm">
                <thead className="bg-gray-50 dark:bg-gray-700">
                    <tr>
                        <th className="px-3 py-2 text-left font-medium text-gray-600 dark:text-gray-300">Key</th>
                        <th className="px-3 py-2 text-left font-medium text-gray-600 dark:text-gray-300">Value</th>
                    </tr>
                </thead>
                <tbody className="divide-y divide-gray-200 dark:divide-gray-600">
                    {entries.map(([key, value]) => (
                        <tr key={key} className="hover:bg-gray-50 dark:hover:bg-gray-700">
                            <td className="px-3 py-2 font-mono text-xs text-gray-800 dark:text-gray-200">{key}</td>
                            <td className="px-3 py-2 text-gray-600 dark:text-gray-400">
                                {typeof value === 'object' ? JSON.stringify(value, null, 2).slice(0, 100) + '...' : String(value)}
                            </td>
                        </tr>
                    ))}
                </tbody>
            </table>
            {Object.entries(data).length > maxRows && (
                <div className="text-xs text-gray-500 dark:text-gray-400 mt-2 text-center">
                    Showing {maxRows} of {Object.entries(data).length} entries
                </div>
            )}
        </div>
    );
};

const Spinner = () => {
    return (
        <div
            className="w-16 h-16 border-8 mt-10 border-blue-500 border-dashed rounded-full animate-spin"
            style={{
                borderColor: "blue transparent blue transparent",
            }}
        />
    );
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

const Message = ({ children }: PropsWithChildren) => (
    <p className="text-lg font-bold text-[#333] dark:text-white text-center">
        {children}
    </p>
);

const ErrorContainer = ({ children }: PropsWithChildren) => (
    <div className="bg-[#fff0f0] text-[#d8000c] border border-[#ffbaba] p-5 my-5 rounded shadow-sm text-left w-[90%] max-w-[600px]">
        {children}
    </div>
);

const ErrorHeading = ({ children }: PropsWithChildren) => (
    <h4 className="mt-0 text-lg font-bold">{children}</h4>
);

// Types for Zarr metadata
interface ZarrMetadata {
    datasetStructure: {
        groups: string[];
        arrays: string[];
    };
    spatialData?: {
        coordinateSystems: Record<string, any>;
        elements: Record<string, any>;
        transformations: any[];
    };
    images?: {
        [key: string]: {
            shape: number[];
            dtype: string;
            chunks: number[];
            dimensions: string[];
            physicalSizeX?: number;
            physicalSizeY?: number;
            physicalSizeZ?: number;
            units?: string;
            channels?: string[];
        };
    };
    tables?: {
        [key: string]: {
            shape: number[];
            obsNames: string[];
            varNames: string[];
            observations: Record<string, any>;
            variables: Record<string, any>;
            embeddings?: string[];
        };
    };
    labels?: {
        [key: string]: {
            shape: number[];
            dtype: string;
            labelValues: number[];
            categories?: string[];
        };
    };
    groupAttributes: Record<string, any>;
}

type MetadataActionType =
    | "SET_URL"
    | "SET_IS_LOADING"
    | "SET_SUCCESS"
    | "SET_ERROR"
    | "SET_METADATA";

const DEFAULT_REDUCER_STATE = {
    url: "",
    isLoading: false,
    success: false,
    error: null as unknown,
    metadata: null as ZarrMetadata | null,
} as const;

type ReducerState = typeof DEFAULT_REDUCER_STATE;

type ReducerPayload<T extends MetadataActionType> = T extends "SET_URL"
    ? string
    : T extends "SET_IS_LOADING"
      ? boolean
      : T extends "SET_SUCCESS"
        ? boolean
        : T extends "SET_ERROR"
          ? { message: string; details?: string }
          : T extends "SET_METADATA"
            ? ZarrMetadata
            : never;

const reducer = <T extends MetadataActionType>(
    state: ReducerState,
    action: { type: T; payload: any },
) => {
    switch (action.type) {
        case "SET_URL":
            return { ...state, url: action.payload };
        case "SET_IS_LOADING":
            return { ...state, isLoading: action.payload };
        case "SET_SUCCESS":
            return { ...state, success: action.payload };
        case "SET_ERROR":
            return { ...state, error: action.payload };
        case "SET_METADATA":
            return { ...state, metadata: action.payload };
        default:
            return state;
    }
};

export interface ZarrMetadataDialogComponentProps {
    onClose: () => void;
    onResize: (width: number, height: number) => void;
}

const ZarrMetadataDialogComponent: React.FC<ZarrMetadataDialogComponentProps> =
    observer(({ onClose, onResize }: ZarrMetadataDialogComponentProps) => {
        const { root } = useProject();
        const [state, dispatch] = useReducer(reducer, DEFAULT_REDUCER_STATE);

        // Fetch zarr metadata from backend API
        const fetchZarrMetadata = async (url: string): Promise<ZarrMetadata> => {
            const apiUrl = `${root}/zarr/metadata?url=${encodeURIComponent(url)}`;
            
            console.log('ZarrMetadataDialog DEBUG:');
            console.log('- root:', root);
            console.log('- constructed apiUrl:', apiUrl);
            console.log('- zarr url parameter:', url);
            
            const response = await fetch(apiUrl, {
                method: 'GET',
                headers: {
                    'Content-Type': 'application/json',
                },
            });

            if (!response.ok) {
                const errorData = await response.json().catch(() => ({}));
                throw new Error(errorData.message || `HTTP error! status: ${response.status}`);
            }

            const data = await response.json();
            
            if (!data.success) {
                throw new Error(data.message || 'Failed to fetch metadata');
            }

            return data.metadata;
        };

        const handleUrlChange = useCallback((event: React.ChangeEvent<HTMLInputElement>) => {
            dispatch({ type: "SET_URL", payload: event.target.value });
        }, []);

        const handleFetchMetadata = async () => {
            if (!state.url.trim()) {
                dispatch({
                    type: "SET_ERROR",
                    payload: {
                        message: "Please enter a valid Zarr dataset URL",
                    },
                });
                return;
            }

            dispatch({ type: "SET_IS_LOADING", payload: true });
            dispatch({ type: "SET_ERROR", payload: null });

            try {
                const zarrUrl = state.url.trim();
                
                // Call backend API to fetch metadata
                const metadata = await fetchZarrMetadata(zarrUrl);

                if (!metadata) {
                    throw new Error("Metadata could not be extracted.");
                }

                dispatch({ type: "SET_METADATA", payload: metadata });
                dispatch({ type: "SET_SUCCESS", payload: true });
                onResize(1200, 800); 

            } catch (error) {
                console.error("Error fetching Zarr metadata:", error);
                dispatch({
                    type: "SET_ERROR",
                    payload: {
                        message: "Failed to fetch metadata from the provided URL",
                        details: error instanceof Error ? error.message : String(error),
                    },
                });
            } finally {
                dispatch({ type: "SET_IS_LOADING", payload: false });
            }
        };

        const handleClose = useCallback(() => {
            onResize(600, 400); // Reset to default size
            onClose();
        }, [onClose, onResize]);

        const resetForm = () => {
            dispatch({ type: "SET_URL", payload: "" });
            dispatch({ type: "SET_ERROR", payload: null });
            dispatch({ type: "SET_SUCCESS", payload: false });
            dispatch({ type: "SET_METADATA", payload: null });
            onResize(600, 400);
        };

        // Set initial dialog size
        useEffect(() => {
            onResize(600, 400);
        }, [onResize]);

        const renderMetadataContent = () => {
            if (!state.metadata) return null;

            const { metadata } = state;
            const attrs = metadata.groupAttributes;
            
            // Extract data from different groups
            const cellsAttrs = attrs.cells || {};
            const transcriptsAttrs = attrs.transcripts || {};
            const analysisAttrs = attrs.analysis || {};
            const rootAttrs = { ...attrs };
            delete rootAttrs.cells;
            delete rootAttrs.transcripts;
            delete rootAttrs.analysis;

            return (
                <div className="w-full max-h-[600px] overflow-y-auto">
                    {/* Dataset Overview */}
                    <MetadataContainer>
                        <SectionHeader icon={<InfoIcon />}>Dataset Overview</SectionHeader>
                        <MetadataGrid>
                            <MetadataCard title="Dataset URL">
                                <span className="font-mono text-xs break-all">{state.url}</span>
                            </MetadataCard>
                            <MetadataCard title="Zarr Structure">
                                <div className="space-y-1 text-xs">
                                    <div><strong>Groups:</strong> {metadata.datasetStructure.groups.length}</div>
                                    <div><strong>Arrays:</strong> {metadata.datasetStructure.arrays.length}</div>
                                    <div><strong>Format:</strong> Zarr v2</div>
                                </div>
                            </MetadataCard>
                            <MetadataCard title="Data Tables">
                                <div className="space-y-1 text-xs">
                                    {Object.keys(metadata.tables || {}).length > 0 ? (
                                        Object.entries(metadata.tables || {}).map(([name, info]: [string, any]) => (
                                            <div key={name}>
                                                <strong>{name}:</strong> {info.shape?.join(' × ') || 'Unknown shape'}
                                            </div>
                                        ))
                                    ) : (
                                        <div className="text-gray-500">No table data found</div>
                                    )}
                                </div>
                            </MetadataCard>
                        </MetadataGrid>
                        
                        {/* Discovered Groups */}
                        {metadata.datasetStructure.groups.length > 0 && (
                            <div className="mt-4">
                                <h5 className="font-medium text-sm mb-2">Discovered Zarr Groups</h5>
                                <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-4 gap-2">
                                    {metadata.datasetStructure.groups.map((group: string) => (
                                        <div key={group} className="bg-blue-50 dark:bg-blue-900 px-3 py-2 rounded text-xs font-mono text-center">
                                            {group}
                                        </div>
                                    ))}
                                </div>
                            </div>
                        )}

                        {/* Discovered Arrays */}
                        {metadata.datasetStructure.arrays.length > 0 && (
                            <div className="mt-4">
                                <h5 className="font-medium text-sm mb-2">Discovered Zarr Arrays</h5>
                                <div className="grid grid-cols-1 md:grid-cols-2 gap-2">
                                    {metadata.datasetStructure.arrays.map((array: string) => (
                                        <div key={array} className="bg-green-50 dark:bg-green-900 px-3 py-2 rounded">
                                            <div className="font-mono text-xs font-semibold">{array}</div>
                                            {metadata.tables?.[array.split('/').pop() || ''] && (
                                                <div className="text-xs text-gray-600 dark:text-gray-400 mt-1">
                                                    Shape: {metadata.tables[array.split('/').pop() || ''].shape?.join(' × ')} | 
                                                    Type: {metadata.tables[array.split('/').pop() || ''].dtype}
                                                </div>
                                            )}
                                        </div>
                                    ))}
                                </div>
                            </div>
                        )}
                    </MetadataContainer>

                    {/* Cell Data Section */}
                    {Object.keys(cellsAttrs).length > 0 && (
                        <MetadataContainer>
                            <SectionHeader icon={<TableIcon />}>Cell Data (cells.zarr)</SectionHeader>
                            <MetadataGrid>
                                <MetadataCard title="Cell Arrays">
                                    <div className="space-y-1 text-xs">
                                        {metadata.datasetStructure.arrays
                                            .filter((arr: string) => arr.startsWith('cells.zarr/'))
                                            .map((arr: string) => (
                                                <div key={arr} className="font-mono">{arr.replace('cells.zarr/', '')}</div>
                                            ))
                                        }
                                    </div>
                                </MetadataCard>
                                <MetadataCard title="Cell Attributes">
                                    <div className="space-y-1 text-xs">
                                        {Object.entries(cellsAttrs).slice(0, 5).map(([key, value]) => (
                                            <div key={key}>
                                                <strong>{key}:</strong> {typeof value === 'object' ? JSON.stringify(value).slice(0, 50) + '...' : String(value)}
                                            </div>
                                        ))}
                                        {Object.keys(cellsAttrs).length > 5 && (
                                            <div className="text-gray-500">...and {Object.keys(cellsAttrs).length - 5} more</div>
                                        )}
                                    </div>
                                </MetadataCard>
                            </MetadataGrid>
                        </MetadataContainer>
                    )}

                    {/* Transcript Data Section */}
                    {Object.keys(transcriptsAttrs).length > 0 && (
                        <MetadataContainer>
                            <SectionHeader icon={<TableIcon />}>Transcript Data (transcripts.zarr)</SectionHeader>
                            <MetadataGrid>
                                <MetadataCard title="RNA Counts">
                                    <div className="space-y-1 text-xs">
                                        <div><strong>Total RNAs:</strong> {transcriptsAttrs.number_rnas?.toLocaleString() || 'Unknown'}</div>
                                        <div><strong>Total Genes:</strong> {transcriptsAttrs.number_genes?.toLocaleString() || 'Unknown'}</div>
                                        <div><strong>Codewords:</strong> {transcriptsAttrs.codeword_count?.toLocaleString() || 'Unknown'}</div>
                                    </div>
                                </MetadataCard>
                                <MetadataCard title="Gene Categories">
                                    <div className="space-y-1 text-xs">
                                        {transcriptsAttrs.codeword_gene_names && (
                                            <>
                                                <div><strong>Target Genes:</strong> {transcriptsAttrs.codeword_gene_names.filter((name: string) => !name.startsWith('NegControl') && !name.startsWith('UnassignedCodeword')).length}</div>
                                                <div><strong>Negative Controls:</strong> {transcriptsAttrs.codeword_gene_names.filter((name: string) => name.startsWith('NegControl')).length}</div>
                                                <div><strong>Unassigned:</strong> {transcriptsAttrs.codeword_gene_names.filter((name: string) => name.startsWith('UnassignedCodeword')).length}</div>
                                            </>
                                        )}
                                    </div>
                                </MetadataCard>
                                <MetadataCard title="Transcript Arrays">
                                    <div className="space-y-1 text-xs">
                                        {metadata.datasetStructure.arrays
                                            .filter((arr: string) => arr.startsWith('transcripts.zarr/'))
                                            .slice(0, 8)
                                            .map((arr: string) => (
                                                <div key={arr} className="font-mono">{arr.replace('transcripts.zarr/', '')}</div>
                                            ))
                                        }
                                        {metadata.datasetStructure.arrays.filter((arr: string) => arr.startsWith('transcripts.zarr/')).length > 8 && (
                                            <div className="text-gray-500">...and {metadata.datasetStructure.arrays.filter((arr: string) => arr.startsWith('transcripts.zarr/')).length - 8} more</div>
                                        )}
                                    </div>
                                </MetadataCard>
                            </MetadataGrid>
                            
                            {/* Sample Genes */}
                            {transcriptsAttrs.gene_names && (
                                <div className="mt-4">
                                    <h5 className="font-medium text-sm mb-2">Sample Target Genes</h5>
                                    <div className="grid grid-cols-3 md:grid-cols-6 gap-2 text-xs">
                                        {transcriptsAttrs.gene_names.filter((name: string) => !name.startsWith('NegControl') && !name.startsWith('UnassignedCodeword')).slice(0, 24).map((gene: string) => (
                                            <div key={gene} className="font-mono bg-gray-100 dark:bg-gray-700 px-2 py-1 rounded text-center">{gene}</div>
                                        ))}
                                    </div>
                                    {transcriptsAttrs.gene_names.length > 24 && (
                                        <div className="text-xs text-gray-500 mt-2">
                                            Showing 24 of {transcriptsAttrs.gene_names.filter((name: string) => !name.startsWith('NegControl') && !name.startsWith('UnassignedCodeword')).length} target genes
                                        </div>
                                    )}
                                </div>
                            )}

                            {/* Fields of View */}
                            {transcriptsAttrs.fov_names && (
                                <div className="mt-4">
                                    <h5 className="font-medium text-sm mb-2">Fields of View ({transcriptsAttrs.fov_names.length} total)</h5>
                                    <div className="grid grid-cols-4 md:grid-cols-8 gap-1 text-xs max-h-32 overflow-y-auto">
                                        {transcriptsAttrs.fov_names.slice(0, 32).map((fov: string) => (
                                            <div key={fov} className="font-mono bg-blue-50 dark:bg-blue-900 px-1 py-1 rounded text-center">{fov}</div>
                                        ))}
                                    </div>
                                    {transcriptsAttrs.fov_names.length > 32 && (
                                        <div className="text-xs text-gray-500 mt-2">
                                            Showing first 32 fields of view
                                        </div>
                                    )}
                                </div>
                            )}
                        </MetadataContainer>
                    )}

                    {/* Analysis Data Section */}
                    {Object.keys(analysisAttrs).length > 0 && (
                        <MetadataContainer>
                            <SectionHeader icon={<TableIcon />}>Analysis Data (analysis.zarr)</SectionHeader>
                            <MetadataGrid>
                                <MetadataCard title="Analysis Arrays">
                                    <div className="space-y-1 text-xs">
                                        {metadata.datasetStructure.arrays
                                            .filter((arr: string) => arr.startsWith('analysis.zarr/'))
                                            .map((arr: string) => (
                                                <div key={arr} className="font-mono">{arr.replace('analysis.zarr/', '')}</div>
                                            ))
                                        }
                                    </div>
                                </MetadataCard>
                                <MetadataCard title="Cell Groups">
                                    <div className="space-y-1 text-xs">
                                        {metadata.tables?.cell_groups && (
                                            <>
                                                <div><strong>Shape:</strong> {metadata.tables.cell_groups.shape?.join(' × ') || 'Unknown'}</div>
                                                <div><strong>Type:</strong> {metadata.tables.cell_groups.dtype || 'Unknown'}</div>
                                                <div><strong>Description:</strong> Cell grouping/clustering analysis results</div>
                                            </>
                                        )}
                                    </div>
                                </MetadataCard>
                                <MetadataCard title="Analysis Attributes">
                                    <div className="space-y-1 text-xs">
                                        {Object.entries(analysisAttrs).slice(0, 5).map(([key, value]) => (
                                            <div key={key}>
                                                <strong>{key}:</strong> {typeof value === 'object' ? JSON.stringify(value).slice(0, 50) + '...' : String(value)}
                                            </div>
                                        ))}
                                        {Object.keys(analysisAttrs).length > 5 && (
                                            <div className="text-gray-500">...and {Object.keys(analysisAttrs).length - 5} more</div>
                                        )}
                                    </div>
                                </MetadataCard>
                            </MetadataGrid>
                        </MetadataContainer>
                    )}

                    {/* Additional Root Attributes */}
                    {Object.keys(rootAttrs).length > 0 && (
                        <MetadataContainer>
                            <SectionHeader icon={<FolderIcon />}>Root Attributes</SectionHeader>
                            <div className="space-y-2">
                                <DataTable 
                                    data={Object.fromEntries(Object.entries(rootAttrs).filter(([key, value]) => 
                                        !['codeword_gene_mapping', 'codeword_gene_names', 'gene_names', 'gene_index_map', 'fov_names', 'cell_groups'].includes(key)
                                    ).slice(0, 15))} 
                                    maxRows={15} 
                                />
                                {Object.keys(rootAttrs).length > 15 && (
                                    <div className="text-xs text-gray-500 mt-2">
                                        Showing 15 of {Object.keys(rootAttrs).length} root attributes
                                    </div>
                                )}
                            </div>
                        </MetadataContainer>
                    )}

                    {/* Gene Index Mapping */}
                    {(transcriptsAttrs.gene_index_map || rootAttrs.gene_index_map) && (
                        <MetadataContainer>
                            <SectionHeader icon={<TableIcon />}>Gene Index Mapping</SectionHeader>
                            <div className="space-y-2">
                                <div className="text-sm text-gray-600 dark:text-gray-400 mb-2">
                                    Maps gene names to their indices in the expression matrix
                                </div>
                                <DataTable 
                                    data={Object.fromEntries(
                                        Object.entries(transcriptsAttrs.gene_index_map || rootAttrs.gene_index_map || {})
                                            .filter(([gene]) => !gene.startsWith('NegControl') && !gene.startsWith('UnassignedCodeword'))
                                            .slice(0, 25)
                                    )} 
                                    maxRows={25} 
                                />
                                <div className="text-xs text-gray-500 mt-2">
                                    Showing sample gene mappings (target genes only)
                                </div>
                            </div>
                        </MetadataContainer>
                    )}
                </div>
            );
        };
        return (
            <Container>
                {state.isLoading ? (
                    <StatusContainer>
                        <Message>Fetching Zarr metadata, please wait...</Message>
                        <Spinner />
                    </StatusContainer>
                ) : state.error ? (
                    <>
                        <ErrorContainer>
                            <ErrorHeading>Error</ErrorHeading>
                            <p>{(state.error as any).message}</p>
                            {(state.error as any).details && (
                                <details className="mt-2">
                                    <summary className="cursor-pointer text-sm">Details</summary>
                                    <pre className="text-xs mt-1 p-2 bg-gray-100 dark:bg-gray-800 rounded overflow-auto">
                                        {(state.error as any).details}
                                    </pre>
                                </details>
                            )}
                        </ErrorContainer>
                        <div className="flex justify-center items-center gap-6">
                            <Button marginTop="mt-1" onClick={resetForm}>
                                Try Again
                            </Button>
                            <Button
                                color="red"
                                size="px-6 py-2.5"
                                marginTop="mt-1"
                                onClick={handleClose}
                            >
                                Close
                            </Button>
                        </div>
                    </>
                ) : state.success && state.metadata ? (
                    <>
                        {renderMetadataContent()}
                        <div className="flex justify-center items-center gap-6 mt-6">
                            <Button marginTop="mt-1" onClick={resetForm}>
                                Link Dataset to Project
                            </Button>
                            <Button
                                color="red"
                                size="px-6 py-2.5"
                                marginTop="mt-1"
                                onClick={handleClose}
                            >
                                Close
                            </Button>
                        </div>
                    </>
                ) : (
                    <>
                        <div className="flex flex-col items-center justify-center space-y-6 p-6">
                            <CloudDownloadIcon
                                className="text-blue-500"
                                style={{ fontSize: "4rem" }}
                            />
                            <div className="w-full max-w-lg">
                                <label htmlFor="zarr-dataset-url" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
                                    Xenium Zarr Dataset URL
                                </label>
                                <input
                                    id="zarr-dataset-url"
                                    type="url"
                                    value={state.url}
                                    onChange={handleUrlChange}
                                    className="w-full p-3 border rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 dark:text-gray-300 dark:bg-gray-800 dark:border-gray-600"
                                    placeholder="http://localhost:8000/your-zarr-dataset.zarr"
                                />
                                <p className="text-xs text-gray-500 dark:text-gray-400 mt-1">
                                    Enter the URL of a Zarr dataset (supports HTTP/HTTPS)
                                </p>
                            </div>
                            <div className="flex justify-center items-center gap-6">
                                <Button 
                                    onClick={handleFetchMetadata}
                                    disabled={!state.url.trim()}
                                >
                                    Fetch Metadata
                                </Button>
                                <Button
                                    color="gray"
                                    onClick={handleClose}
                                >
                                    Cancel
                                </Button>
                            </div>
                        </div>
                    </>
                )}
            </Container>
        );
    });

const Wrapper = (props: ZarrMetadataDialogComponentProps) => {
    const [open, setOpen] = useState(true);

    const handleClose = () => {
        setOpen(false);
        props.onClose();
    };

    return (
        <ReusableDialog 
            open={open} 
            handleClose={handleClose} 
            component={<ZarrMetadataDialogComponent {...props} onClose={handleClose} />} 
        />
    );
};

export default Wrapper;