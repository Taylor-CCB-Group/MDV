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

        // Mock function to simulate zarr metadata fetching
        // Replace this with actual zarr-js implementation
        const fetchZarrMetadata = async (url: string): Promise<ZarrMetadata> => {
            // Simulate network delay
            await new Promise(resolve => setTimeout(resolve, 2000));
            
            // Mock data structure based on Xenium dataset
            return {
                datasetStructure: {
                    groups: ["images", "labels", "tables"],
                    arrays: ["images/tissue_image", "labels/nuclei_segmentation", "tables/my_ann_data"]
                },
                spatialData: {
                    coordinateSystems: {
                        "global": {
                            "axes": ["x", "y"],
                            "units": ["micrometer", "micrometer"]
                        }
                    },
                    elements: {
                        "tissue_image": "image",
                        "nuclei_segmentation": "labels",
                        "my_ann_data": "table"
                    },
                    transformations: []
                },
                images: {
                    "tissue_image": {
                        shape: [4096, 4096, 3],
                        dtype: "uint8",
                        chunks: [512, 512, 3],
                        dimensions: ["y", "x", "c"],
                        physicalSizeX: 0.2125,
                        physicalSizeY: 0.2125,
                        units: "micrometer",
                        channels: ["Red", "Green", "Blue"]
                    }
                },
                tables: {
                    "my_ann_data": {
                        shape: [156432, 541],
                        obsNames: ["cell_001", "cell_002", "cell_003"],
                        varNames: ["GAPDH", "ACTB", "MYC"],
                        observations: {
                            "cell_type": "string",
                            "total_counts": "float64",
                            "n_genes_by_counts": "int64"
                        },
                        variables: {
                            "gene_symbol": "string",
                            "highly_variable": "bool",
                            "mean": "float64"
                        },
                        embeddings: ["X_umap", "X_pca"]
                    }
                },
                labels: {
                    "nuclei_segmentation": {
                        shape: [4096, 4096],
                        dtype: "uint32",
                        labelValues: [0, 1, 2, 3],
                        categories: ["background", "nucleus"]
                    }
                },
                groupAttributes: {
                    "created_by": "spatialdata-io",
                    "version": "0.1.0",
                    "experiment_type": "Xenium",
                    "organism": "Homo sapiens"
                }
            };
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
                let transformedUrl = state.url.trim();

                // Match EMBL Zarr URLs and transform them for the local proxy
                if (transformedUrl.startsWith("https://s3.embl.de/spatialdata/")) {
                    const subpath = transformedUrl.replace("https://s3.embl.de/spatialdata/", "");
                    transformedUrl = `${window.location.origin}/zarr_proxy/${subpath}`;
                } else if (transformedUrl.startsWith('/zarr_proxy/')) {
                    transformedUrl = `${window.location.origin}${transformedUrl}`;
                }
                
                // For debugging, you can keep this, but it's not needed for the final logic.
                // await debugExploreZarrDataset(transformedUrl);

                // Call the primary extraction function. It has its own internal fallback.
                const metadata = await extractZarrMetadata(transformedUrl);

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

            return (
                <div className="w-full max-h-[600px] overflow-y-auto">
                    {/* Dataset Overview */}
                    <MetadataContainer>
                        <SectionHeader icon={<InfoIcon />}>Dataset Overview</SectionHeader>
                        <MetadataGrid>
                            <MetadataCard title="Dataset URL">
                                <span className="font-mono text-xs break-all">{state.url}</span>
                            </MetadataCard>
                            <MetadataCard title="Dataset Info">
                                <div className="space-y-1 text-xs">
                                    <div><strong>Name:</strong> {attrs.name || 'Unknown'}</div>
                                    <div><strong>Version:</strong> {attrs.major_version}.{attrs.minor_version}</div>
                                    <div><strong>UUID:</strong> {attrs.dataset_uuid?.slice(0, 8)}...</div>
                                    <div><strong>Format:</strong> {attrs.data_format}</div>
                                </div>
                            </MetadataCard>
                            <MetadataCard title="Spatial Info">
                                <div className="space-y-1 text-xs">
                                    <div><strong>Units:</strong> {attrs.spatial_units}</div>
                                    <div><strong>Coordinate Space:</strong> {attrs.coordinate_space}</div>
                                    <div><strong>FOV Count:</strong> {attrs.fov_names?.length || 0}</div>
                                </div>
                            </MetadataCard>
                        </MetadataGrid>
                    </MetadataContainer>

                    {/* Transcripts Data */}
                    <MetadataContainer>
                        <SectionHeader icon={<TableIcon />}>Transcript Data</SectionHeader>
                        <MetadataGrid>
                            <MetadataCard title="RNA Counts">
                                <div className="space-y-1 text-xs">
                                    <div><strong>Total RNAs:</strong> {attrs.number_rnas?.toLocaleString() || 'Unknown'}</div>
                                    <div><strong>Total Genes:</strong> {attrs.number_genes?.toLocaleString() || 'Unknown'}</div>
                                    <div><strong>Codewords:</strong> {attrs.codeword_count?.toLocaleString() || 'Unknown'}</div>
                                </div>
                            </MetadataCard>
                            <MetadataCard title="Gene Categories">
                                <div className="space-y-1 text-xs">
                                    {attrs.codeword_gene_names && (
                                        <>
                                            <div><strong>Target Genes:</strong> {attrs.codeword_gene_names.filter((name: string) => !name.startsWith('NegControl') && !name.startsWith('UnassignedCodeword')).length}</div>
                                            <div><strong>Negative Controls:</strong> {attrs.codeword_gene_names.filter((name: string) => name.startsWith('NegControl')).length}</div>
                                            <div><strong>Unassigned:</strong> {attrs.codeword_gene_names.filter((name: string) => name.startsWith('UnassignedCodeword')).length}</div>
                                        </>
                                    )}
                                </div>
                            </MetadataCard>
                            <MetadataCard title="Fields of View">
                                <div className="space-y-1 text-xs max-h-20 overflow-y-auto">
                                    {attrs.fov_names?.slice(0, 10).map((fov: string) => (
                                        <div key={fov} className="font-mono text-xs">{fov}</div>
                                    ))}
                                    {attrs.fov_names?.length > 10 && <div className="text-gray-500">...and {attrs.fov_names.length - 10} more</div>}
                                </div>
                            </MetadataCard>
                        </MetadataGrid>
                        
                        {/* Sample Genes */}
                        {attrs.gene_names && (
                            <div className="mt-4">
                                <h5 className="font-medium text-sm mb-2">Sample Target Genes</h5>
                                <div className="grid grid-cols-3 md:grid-cols-6 gap-2 text-xs">
                                    {attrs.gene_names.filter((name: string) => !name.startsWith('NegControl') && !name.startsWith('UnassignedCodeword')).slice(0, 18).map((gene: string) => (
                                        <div key={gene} className="font-mono bg-gray-100 dark:bg-gray-700 px-2 py-1 rounded text-center">{gene}</div>
                                    ))}
                                </div>
                                {attrs.gene_names.length > 18 && (
                                    <div className="text-xs text-gray-500 mt-2">
                                        Showing 18 of {attrs.gene_names.filter((name: string) => !name.startsWith('NegControl') && !name.startsWith('UnassignedCodeword')).length} target genes
                                    </div>
                                )}
                            </div>
                        )}
                    </MetadataContainer>

                    {/* Raw Structure Explorer */}
                    {metadata.rawStructure && (
                        <MetadataContainer>
                            <SectionHeader icon={<FolderIcon />}>Zarr Structure</SectionHeader>
                            <div className="space-y-2">
                                <div className="text-sm">
                                    <strong>Groups:</strong> {metadata.datasetStructure.groups.length} | <strong>Arrays:</strong> {metadata.datasetStructure.arrays.length}
                                </div>
                                {metadata.datasetStructure.groups.length === 0 && metadata.datasetStructure.arrays.length === 0 && (
                                    <div className="text-xs text-gray-600 dark:text-gray-400 bg-gray-100 dark:bg-gray-700 p-2 rounded">
                                        This appears to be a root-level Zarr group with metadata only. The subdirectories (codeword_category, density, gene_category, grids) were not accessible as Zarr arrays or groups.
                                    </div>
                                )}
                                {Object.keys(attrs).length > 0 && (
                                    <div>
                                        <h5 className="font-medium text-sm mb-2">Root Attributes ({Object.keys(attrs).length} total)</h5>
                                        <DataTable data={Object.fromEntries(Object.entries(attrs).filter(([key, value]) => 
                                            !['codeword_gene_mapping', 'codeword_gene_names', 'gene_names', 'gene_index_map', 'fov_names'].includes(key)
                                        ).slice(0, 10))} maxRows={10} />
                                    </div>
                                )}
                            </div>
                        </MetadataContainer>
                    )}

                    {/* Gene Index Mapping */}
                    {attrs.gene_index_map && (
                        <MetadataContainer>
                            <SectionHeader icon={<TableIcon />}>Gene Index Mapping</SectionHeader>
                            <div className="space-y-2">
                                <div className="text-sm text-gray-600 dark:text-gray-400 mb-2">
                                    Maps gene names to their indices in the expression matrix
                                </div>
                                <DataTable 
                                    data={Object.fromEntries(
                                        Object.entries(attrs.gene_index_map)
                                            .filter(([gene]) => !gene.startsWith('NegControl') && !gene.startsWith('UnassignedCodeword'))
                                            .slice(0, 20)
                                    )} 
                                    maxRows={20} 
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
                                    placeholder="https://s3.embl.de/spatialdata/spatialdata-sandbox/xenium_rep1_io.zarr"
                                />
                                <p className="text-xs text-gray-500 dark:text-gray-400 mt-1">
                                    Enter the URL of a remote Zarr dataset in SpatialData format
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