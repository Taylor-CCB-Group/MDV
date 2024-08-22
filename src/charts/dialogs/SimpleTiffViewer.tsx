import React, { useEffect, useMemo } from 'react';
import { observer } from "mobx-react-lite";
import { DetailView, getDefaultInitialViewState } from '@hms-dbmi/viv';
import { ColorPaletteExtension } from '@hms-dbmi/viv';
import { useLoader, useChannelsStore, useViewerStore, type VivContextType, useViewerStoreApi } from '../../react/components/avivatorish/state';
import { shallow } from 'zustand/shallow';
import MDVivViewer from '../../react/components/avivatorish/MDVivViewer';
import VivContrastExtension from "@/webgl/VivContrastExtension";

interface SimpleTiffViewerProps {
    width: number;
    height: number;
  }
  
  const SimpleTiffViewer: React.FC<SimpleTiffViewerProps> = observer(({ width, height }) => {
    const loader = useLoader();
    const { colors, contrastLimits, channelsVisible, selections, brightness, contrast } = useChannelsStore(
      ({ colors, contrastLimits, channelsVisible, selections, brightness, contrast }) => ({
        colors, contrastLimits, channelsVisible, selections, brightness, contrast
      }),
      shallow
    );
    const viewState = useViewerStore(store => store.viewState);
  
    useEffect(() => {
      console.log('Loader:', loader);
      console.log('ViewState:', viewState);
    }, [loader, viewState]);
    
    const extensions = useMemo(() => [new ColorPaletteExtension()], []);
    
    const layerConfig = useMemo(() => {
        const config = {
          loader,
          selections,
          contrastLimits,
          colors,
          channelsVisible,
          brightness,
          contrast,
          extensions
        };
        console.log('Layer Config:', config);
        return config;
      }, [loader, selections, contrastLimits, colors, channelsVisible, brightness, contrast, extensions]);
  
  
    const viewerStore = useViewerStoreApi();
    useEffect(() => {
      if (loader && !viewState) {
        const initialViewState = getDefaultInitialViewState(loader, { width, height });
        viewerStore.setState({ viewState: initialViewState });
      }
    }, [loader, viewState, width, height, viewerStore]);
  
    const detailView = useMemo(() => new DetailView({
      id: 'detail-view',
      snapScaleBar: true,
      width,
      height
    }), [width, height]);
    
  
    if (!loader) {
      return <div>No loader available</div>;
    }
  
    if (!viewState) {
      return <div>No view state available</div>;
    }
  
    return (
        <MDVivViewer
          views={[detailView]}
          layerProps={[layerConfig]}
          viewStates={[{ ...viewState, id: 'detail-view' }]}
          onViewStateChange={e => {
            console.log('View State Change:', e.viewState);
            viewerStore.setState({ viewState: { ...e.viewState, id: 'detail-view' } });
          }}
          deckProps={{
            glOptions: {
              preserveDrawingBuffer: true,
            },
            style: {
              height: '200px',
              width: '200px',
              position: 'relative',
            },
            onError: (error) => {
              console.error('MDVivViewer Error:', error);
            },
          }}
        />
      );
    });
  
  export default SimpleTiffViewer;