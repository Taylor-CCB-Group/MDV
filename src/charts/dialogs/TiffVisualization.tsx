import React, { useEffect, useMemo, useState } from 'react';
import { observer } from 'mobx-react-lite';
import { useViewerStoreApi, useViewerStore, } from '../../react/components/avivatorish/state';
import type { OME_TIFF } from '../../react/components/avivatorish/state';
import SimpleTiffViewer from './SimpleTiffViewer';
import { useImage } from '@/react/components/avivatorish/hooks';
import TiffChartWrapper from './TiffChartWrapper';

const TiffVisualization = observer(({ metadata, file }: { metadata: OME_TIFF['metadata']; file: File; }) => {
  const [fileUrl, setFileUrl] = useState<string | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [isFullscreen, setIsFullscreen] = useState(false);
  
  const viewerStore = useViewerStoreApi();

  useEffect(() => {
    if (file) {
      const url = URL.createObjectURL(file);
      setFileUrl(url);

      return () => {
        URL.revokeObjectURL(url);
      };
    }
  }, [file]);

  useEffect(() => {
    if (fileUrl) {
      setIsLoading(true);
      const source: { urlOrFile: string; description: string } = { 
        urlOrFile: fileUrl, 
        description: file.name 
      };
      viewerStore.setState({ source });
      
      setTimeout(() => {
        setIsLoading(false);
      }, 1000);
    }
  }, [fileUrl, file.name, viewerStore]);

  const source = useViewerStore(store => store.source);
  useImage(source);

  const handleFullscreenChange = (fullscreenState: boolean) => {
    setIsFullscreen(fullscreenState);
  };

  if (isLoading) {
    return <div>Loading TIFF data...</div>;
  }

  const config = {
    title: file.name,
    legend: "TIFF Visualization",
  };

  return (
    <TiffChartWrapper 
      metadata={metadata} 
      file={file} 
      config={config} 
      onFullscreenChange={handleFullscreenChange}
    >
      <div className="flex items-center justify-center">
        <SimpleTiffViewer width={200} height={200} isFullscreen={isFullscreen} />
      </div>
    </TiffChartWrapper>
  );
});

export default TiffVisualization;