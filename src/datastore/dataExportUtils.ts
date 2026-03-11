import type { ColumnName } from "@/charts/charts";
import type DataStore from "@/datastore/DataStore";
import { DataModel } from "@/table/DataModel";

/**
 * Create a ReadableStream that will provide data incrementally from the rows and columns
 * selected in the dataModel.
 * @param includeIndex - prepend index column if true, else ignore it
 */
export async function getExportCsvStream(
    dataModel: DataModel,
    delimiter = "\t",
    newline = "\n",
    includeIndex = true
) {
    const columns: ColumnName[] = dataModel.columns;
    const dataStore = dataModel.dataStore;
    dataModel.updateModel(); // Side-effect: updates dataModel.data
    const indexes = dataModel.data;
    const len = indexes.length;
    const encoder = new TextEncoder();

    // Create a ReadableStream that will provide the data incrementally
    const stream = new ReadableStream({
        async start(controller) {
            try {
                // Write the headers first
                const header = (includeIndex ? ["index", ...columns] : columns).join(delimiter) + newline;
                controller.enqueue(encoder.encode(header)); // Add header to stream

                // Write each row of data incrementally
                for (let i = 0; i < len; i++) {
                    const index = indexes[i];
                    const o = dataStore.getRowAsObject(index, columns);
                    const rowValues = columns.map((x) => o[x].toString());
                    // prepend index column if includeIndex is true
                    const line = (includeIndex ? [i.toString(), ...rowValues] : rowValues).join(delimiter) + newline;

                    // Enqueue each encoded line to the stream
                    controller.enqueue(encoder.encode(line));

                    // Optionally add a short delay (yield control) to avoid blocking the main thread
                    // todo: also update a progress bar here
                    if (i % 1e5 === 0) {
                        await new Promise(resolve => setTimeout(resolve, 0));
                    }
                }

                // Close the stream when finished
                controller.close();
                alert("Export complete");
            } catch (err) {
                console.error('Stream writing failed:', err);
                controller.error(err); // Signal an error in the stream
            }
        }
    });

    return stream; // Return the ReadableStream
}

export type TableExportOptions = {
    includeIndex?: boolean;
    delimiter?: string;
    newline?: string;
};

/**
 * Export table data as a blob. Creates a temporary DataModel and uses includeIndex to include index column.
 */
export async function getTableExportBlob(
    dataStore: DataStore,
    columns: string[],
    options: TableExportOptions = {}
): Promise<Blob> {
    const { includeIndex = true, delimiter = "\t", newline = "\n" } = options;
    const dataModel = new DataModel(dataStore, { autoupdate: false });
    dataModel.setColumns(columns);
    const stream = await getExportCsvStream(dataModel, delimiter, newline, includeIndex);
    const response = new Response(stream);
    return response.blob();
}
