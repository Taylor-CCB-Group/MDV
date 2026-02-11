import type React from 'react';
import { useEffect, useMemo } from 'react';
import { observer } from "mobx-react-lite";
import { DetailView, getDefaultInitialViewState } from '@hms-dbmi/viv';
import { ColorPaletteExtension } from '@hms-dbmi/viv';
import { useLoader, useChannelsStore, useViewerStore, useViewerStoreApi } from '../../react/components/avivatorish/state';
import { shallow } from 'zustand/shallow';
import MDVivViewer from '../../react/components/avivatorish/MDVivViewer';
import { useOuterContainer } from '@/react/screen_state';

export interface SimpleTiffViewerProps {
  width: number;
  height: number;
  isFullscreen: boolean;
}

const SimpleTiffViewer: React.FC<SimpleTiffViewerProps> = observer(({ width, height, isFullscreen }: SimpleTiffViewerProps) => {
  const loader = useLoader();
  const { colors, contrastLimits, channelsVisible, selections, brightness, contrast } = useChannelsStore(
    ({ colors, contrastLimits, channelsVisible, selections, brightness, contrast }) => ({
      colors, contrastLimits, channelsVisible, selections, brightness, contrast
    }),
    shallow
  );
  const viewState = useViewerStore(store => store.viewState);
  const outerContainer = useOuterContainer();

  const viewerWidth = isFullscreen ? outerContainer.clientWidth : width;
  const viewerHeight = isFullscreen ? outerContainer.clientHeight : height;

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
      const initialViewState = getDefaultInitialViewState(loader, { width: viewerWidth, height: viewerHeight });
      viewerStore.setState({ viewState: initialViewState });
    }
  }, [loader, viewState, viewerWidth, viewerHeight, viewerStore]);

  const detailView = useMemo(() => new DetailView({
    id: 'detail-view',
    snapScaleBar: true,
    width: viewerWidth,
    height: viewerHeight
  }), [viewerWidth, viewerHeight]);

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
      onViewStateChange={(e: any) => {
        console.log('View State Change:', e.viewState);
        viewerStore.setState({ viewState: { ...e.viewState, id: 'detail-view' } });
      }}
      deckProps={{
        style: {
          height: `${viewerHeight}px`,
          width: `${viewerWidth}px`,
          position: 'relative',
        },
        onError: (error: any) => {
          console.error('MDVivViewer Error:', error);
        },
      }}
    />
  );
});

export default SimpleTiffViewer;