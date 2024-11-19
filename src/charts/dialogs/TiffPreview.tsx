import type React from 'react';
import XMLViewer from 'react-xml-viewer';

// if we could be confident that this was right in general, it would be useful to have at the 
// viv/avivator(ish) level - however, I suspect that there exist other forms of metadata...
// we should attempt to document and characterize these; this is certainly a useful part of that.
// !experimentally using in the context of `useMetadata()` in avivatorish, let's see how it goes.
export type TiffPreviewProps = {
  metadata: {
    ID: string;
    Name: string;
    AquisitionDate: string;
    Description: Record<string, unknown>;
    Pixels: {
      Channels: Array<{
        ID: string;
        SamplesPerPixel: number;
        Name: string;
      }>;
      ID: string;
      DimensionOrder: string;
      Type: string;
      SizeT: number;
      SizeC: number;
      SizeZ: number;
      SizeY: number;
      SizeX: number;
      PhysicalSizeX: number;
      PhysicalSizeY: number;
      SignificantBits: number;
      PhysicalSizeXUnit: string;
      PhysicalSizeYUnit: string;
      PhysicalSizeZUnit: string;
      BigEndian: boolean;
      Interleaved: boolean;
      TiffData: Array<{
        IFD: number;
        PlaneCount: number;
      }>;
    };
  };
}

const jsonToXml = (obj: any, rootName = 'root'): string => {
  const toXml = (obj: any, name: string): string => {
    if (typeof obj === 'object' && obj !== null) {
      if (Array.isArray(obj)) {
        return obj.map(item => toXml(item, name)).join('');
      }
      return `<${name}>${Object.keys(obj).map(key => toXml(obj[key], key)).join('')}</${name}>`;
    }
    return `<${name}>${obj}</${name}>`;
  };

  return `<?xml version="1.0" encoding="UTF-8"?>${toXml(obj, rootName)}`;
};
export const TiffPreview: React.FC<TiffPreviewProps> = ({ metadata }) => {
  const xmlData = jsonToXml(metadata, 'TiffMetadata');

  return (
    <div className="border border-gray-300 rounded-lg shadow-md dark:border-gray-700 w-full mx-10">
      <div className="bg-gray-100 shadow-md rounded-lg p-4 max-w-2xl w-full ">
        <h1 className="text-2xl font-bold mb-4">TIFF Metadata Viewer</h1>
        <div className="bg-white max-h-96 overflow-auto border border-gray-300 rounded-lg p-4">
          <XMLViewer xml={xmlData} initalCollapsedDepth={2} collapsible={true}  />
        </div>
      </div>
    </div>
  );
};
