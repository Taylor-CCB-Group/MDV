interface ColumnPreviewProps {
  columnNames: string[];
  columnTypes: string[];
  secondRowValues: string[];
}

export const ColumnPreview: React.FC<ColumnPreviewProps> = ({ columnNames, columnTypes, secondRowValues }) => {
  return (
    <div className="max-h-80 overflow-y-auto border border-gray-300 rounded-lg shadow-md dark:border-gray-700">
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
          {columnNames.map((columnName, index) => (
            <tr key={columnName}>
              <td className="px-6 py-4 whitespace-nowrap text-sm font-medium text-gray-900 dark:text-gray-100">
                {columnName}
              </td>
              <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-500 dark:text-gray-400">
                {columnTypes[index]}
              </td>
              <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-500 dark:text-gray-400">
                {secondRowValues[index]}
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
};
