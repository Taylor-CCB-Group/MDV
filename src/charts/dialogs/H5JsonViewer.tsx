import { useState } from 'react';
import DebugJsonDialogComponent from '@/react/components/DebugJsonDialogComponent';

export interface JsonViewerProps {
  data: Record<string, any>;
  title?: string;
  initiallyExpanded?: boolean;
  maxHeight?: string;
}

const H5JsonViewer = ({ 
  data, 
  title, 
  initiallyExpanded = false,
  maxHeight = "max-h-96"
}: JsonViewerProps) => {
  const [isExpanded, setIsExpanded] = useState(initiallyExpanded);

  return (
    <div className="bg-white dark:bg-gray-800 rounded-lg shadow-sm border border-gray-200 dark:border-gray-700">
      {title && (
        <div 
          className="flex items-center p-3 cursor-pointer border-b border-gray-200 dark:border-gray-700"
          onClick={() => setIsExpanded(!isExpanded)}
        >
          <h3 className="text-sm font-medium">{title}</h3>
        </div>
      )}
      {isExpanded && (
        <div className={`overflow-auto ${maxHeight}`}>
          <DebugJsonDialogComponent json={data} header={undefined} />
        </div>
      )}
    </div>
  );
};

export default H5JsonViewer;