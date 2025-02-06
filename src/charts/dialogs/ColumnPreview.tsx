import type React from "react";
import { useState } from "react";

export interface ColumnPreviewProps {
    columnNames: string[];
    columnTypes: string[];
    secondRowValues: string[];
}

export const ColumnPreview: React.FC<ColumnPreviewProps> = ({
    columnNames,
    columnTypes,
    secondRowValues,
}) => {
    const [searchQuery, setSearchQuery] = useState<string>("");

    // Filter the columns based on the search query
    const filteredColumns = columnNames
        .map((name, index) => ({
            name,
            type: columnTypes[index],
            value: secondRowValues[index],
        }))
        .filter((column) =>
            column.name.toLowerCase().includes(searchQuery.toLowerCase()),
        );

    return (
        <div className="border border-gray-300 rounded-lg shadow-md dark:border-gray-700 w-full mx-10">
            <div className="p-4 bg-gray-100 dark:bg-gray-900">
                <input
                    type="text"
                    className="w-full px-3 py-2 border border-gray-300 rounded dark:border-gray-700 dark:bg-gray-800 dark:text-white"
                    placeholder="Search column names..."
                    value={searchQuery}
                    onChange={(e) => setSearchQuery(e.target.value)}
                />
            </div>
            <div className="overflow-y-auto h-80 w-full">
                <table className="min-w-full divide-y divide-gray-200 dark:divide-gray-700">
                    <thead className="bg-gray-50 dark:bg-gray-800 sticky top-0 z-10">
                        <tr>
                            <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider dark:text-gray-300">
                                Column Name
                            </th>
                            <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider dark:text-gray-300">
                                Column Type
                            </th>
                            <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider dark:text-gray-300">
                                Second Row Value
                            </th>
                        </tr>
                    </thead>
                    <tbody className="bg-white dark:bg-gray-900 divide-y divide-gray-200 dark:divide-gray-700">
                        {filteredColumns.length > 0 ? (
                            filteredColumns.map((column, index) => (
                                <tr key={index}>
                                    <td className="px-6 py-4 whitespace-nowrap text-sm font-medium text-gray-900 dark:text-gray-100 select-text">
                                        {column.name}
                                    </td>
                                    <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-500 dark:text-gray-400 select-text">
                                        {column.type}
                                    </td>
                                    <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-500 dark:text-gray-400 select-text">
                                        {column.value}
                                    </td>
                                </tr>
                            ))
                        ) : (
                            <tr>
                                <td
                                    colSpan={3}
                                    className="px-6 py-4 text-center text-sm text-gray-500 dark:text-gray-400"
                                >
                                    No columns match your search.
                                </td>
                            </tr>
                        )}
                    </tbody>
                </table>
            </div>
        </div>
    );
};
