import type React from 'react';
import { useState, useMemo } from 'react';

export interface TiffMetadataTableProps {
  metadata: Record<string, any>;
}

export type FlattenedMetadata = {
  element: string;
  attribute: string;
  tag: string;
  value: string;
};

export const flattenObject = (obj: Record<string, any>, prefix: string[] = []): FlattenedMetadata[] => {
  return Object.keys(obj).reduce((acc: FlattenedMetadata[], key) => {
    if (key === 'format') return acc; // Ignore the 'format' key

    const newPrefix = [...prefix, key];
    if (typeof obj[key] === 'object' && obj[key] !== null) {
      if (Array.isArray(obj[key]) && obj[key].every((item: any) => typeof item === 'object')) {
        // Special handling for arrays of objects (like Pixels.Channels)
        obj[key].forEach((item: Record<string, any>, index: number) => {
          Object.keys(item).forEach(itemKey => {
            if (itemKey !== 'format') { // Also ignore 'format' within nested objects
              acc.push({
                element: newPrefix[0] || '',
                attribute: newPrefix[1] || '',
                tag: `${itemKey} (${index})`,
                value: String(item[itemKey] || '')
              });
            }
          });
        });
      } else {
        acc.push(...flattenObject(obj[key], newPrefix));
      }
    } else {
      const [element, attribute] = newPrefix;
      acc.push({
        element: element || '',
        attribute: attribute || '',
        tag: newPrefix.slice(2).join('.') || '',
        value: String(obj[key] || '')
      });
    }
    return acc;
  }, []);
};

export const TiffMetadataTable: React.FC<TiffMetadataTableProps> = ({ metadata }) => {
  const [searchQuery, setSearchQuery] = useState<string>('');

  const flattenedMetadata = useMemo(() => flattenObject(metadata), [metadata]);

  const filteredMetadata = flattenedMetadata.filter(item =>
    (item.element || '').toLowerCase().includes(searchQuery.toLowerCase()) ||
    (item.attribute || '').toLowerCase().includes(searchQuery.toLowerCase()) ||
    (item.tag || '').toLowerCase().includes(searchQuery.toLowerCase()) ||
    (item.value || '').toLowerCase().includes(searchQuery.toLowerCase())
  );

  return (
    <div className="border border-gray-300 rounded-lg shadow-md dark:border-gray-700 w-full mx-10">
      <div className="p-4 bg-gray-100 dark:bg-gray-900">
      <h1 className="text-2xl font-bold mb-4">TIFF Metadata Viewer</h1>
        <input
          type="text"
          className="w-full px-3 py-2 border border-gray-300 rounded dark:border-gray-700 dark:bg-gray-800 dark:text-white"
          placeholder="Search metadata..."
          value={searchQuery}
          onChange={(e) => setSearchQuery(e.target.value)}
        />
      </div>
      <div className="overflow-y-auto h-80 w-full">
        <table className="min-w-full divide-y divide-gray-200 dark:divide-gray-700">
          <thead className="bg-gray-50 dark:bg-gray-800 sticky top-0 z-10">
            <tr>
              <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider dark:text-gray-300 w-1/4">
                Element
              </th>
              <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider dark:text-gray-300 w-1/4">
                Attribute
              </th>
              <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider dark:text-gray-300 w-1/4">
                Tag
              </th>
              <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider dark:text-gray-300 w-1/4">
                Value
              </th>
            </tr>
          </thead>
          <tbody className="bg-white dark:bg-gray-900 divide-y divide-gray-200 dark:divide-gray-700">
            {filteredMetadata.length > 0 ? (
              filteredMetadata.map((item, index) => (
                <tr key={index}>
                  <td className="px-6 py-4 whitespace-nowrap text-sm font-medium text-gray-900 dark:text-gray-100 select-text w-1/4">
                    {item.element || 'N/A'}
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-500 dark:text-gray-400 select-text w-1/4">
                    {item.attribute || 'N/A'}
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-500 dark:text-gray-400 select-text w-1/4">
                    {item.tag || 'N/A'}
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-500 dark:text-gray-400 select-text w-1/4">
                    {item.value || 'N/A'}
                  </td>
                </tr>
              ))
            ) : (
              <tr>
                <td colSpan={4} className="px-6 py-4 text-center text-sm text-gray-500 dark:text-gray-400">
                  No metadata matches your search.
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>
    </div>
  );
};
