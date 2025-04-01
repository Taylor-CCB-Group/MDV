import type { ColumnName } from "@/charts/charts";
import type { DataModel } from "@/table/DataModel";

/**
 * Create a ReadableStream that will provide data incrementally from the rows and columns
 * selected in the dataModel.
 */
export async function getExportCsvStream(
    dataModel: DataModel,
    delimiter = "\t", newline = "\n"
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
                const header = ["index", ...columns].join(delimiter) + newline;
                controller.enqueue(encoder.encode(header)); // Add header to stream

                // Write each row of data incrementally
                for (let i = 0; i < len; i++) {
                    const index = indexes[i];
                    const o = dataStore.getRowAsObject(index, columns);
                    const line = [i.toString()].concat(columns.map((x) => o[x].toString())).join(delimiter) + newline;

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
