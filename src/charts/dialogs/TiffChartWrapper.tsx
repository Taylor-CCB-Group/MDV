import type React from 'react';
import { useState, useCallback, type PropsWithChildren, useRef, useEffect } from 'react'
import { observer } from 'mobx-react-lite';

export interface TiffChartWrapperProps {
  metadata: any;
  file: File;
  config: any;
  onFullscreenChange: (isFullscreen: boolean) => void;
}

const TiffChartWrapper: React.FC<PropsWithChildren<TiffChartWrapperProps>> = observer(({ metadata, file, config, children, onFullscreenChange }
  : PropsWithChildren<TiffChartWrapperProps>) => {
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [oldSize, setOldSize] = useState<[number, number] | null>(null);
  const contentDivRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const handleFullscreenChange = () => {
      const newFullscreenState = !!document.fullscreenElement;
      setIsFullscreen(newFullscreenState);
      onFullscreenChange(newFullscreenState);
      
      if (!newFullscreenState && oldSize) {
        setSize(oldSize[0], oldSize[1]);
        setOldSize(null);
      }
    };

    document.addEventListener('fullscreenchange', handleFullscreenChange);

    return () => {
      document.removeEventListener('fullscreenchange', handleFullscreenChange);
    };
  }, [oldSize, onFullscreenChange]);

  const setSize = (width: number, height: number) => {
    if (contentDivRef.current) {
      contentDivRef.current.style.width = `${width}px`;
      contentDivRef.current.style.height = `${height}px`;
    }
  };

  const toggleFullscreen = async () => {
    if (!isFullscreen) {
      const currentWidth = contentDivRef.current?.offsetWidth || 0;
      const currentHeight = contentDivRef.current?.offsetHeight || 0;
      setOldSize([currentWidth, currentHeight]);
      await contentDivRef.current?.requestFullscreen();
    } else {
      await document.exitFullscreen();
    }
  };

  const addMenuIcon = useCallback((icon: string, tooltip: string, config: any = {}) => {
    return (
      <span
        aria-label={tooltip}
        data-microtip-color="red"
        role="tooltip"
        data-microtip-size={config.size || "small"}
        data-microtip-position={config.position || "bottom-left"}
        className="mx-px cursor-pointer"
        onClick={config.func}
      >
        <i className={`ciview-chart-icon ${icon}`} />
      </span>
    );
  }, []);

  const chartPanelStyle: React.CSSProperties = {
    boxShadow: 'none',
    overflow: 'hidden',
    boxSizing: 'border-box',
    border: '4px solid var(--background_color)',
    width: '100%',
    height: '100%',
  };

  return (
    <div className="w-full h-full" style={chartPanelStyle}>
      <div className="ciview-chart-menuspace flex">
        {addMenuIcon("fas fa-expand", isFullscreen ? "exit fullscreen" : "fullscreen", {
          func: toggleFullscreen
        })}
      </div>
      <div ref={contentDivRef} className="ciview-chart-content flex-grow overflow-auto">
        {children}
      </div>
    </div>
  );
});

export default TiffChartWrapper;
