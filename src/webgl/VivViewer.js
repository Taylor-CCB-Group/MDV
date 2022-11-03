import {
  loadOmeTiff,
  DetailView,
  VolumeView,
  ColorPaletteExtension,
  ColorPalette3DExtensions,
  DETAIL_VIEW_ID,
} from '@hms-dbmi/viv';

import {hexToRGB} from "../datastore/DataStore.js";
import {Deck} from '@deck.gl/core';
import { getMultiSelectionStats, getDefaultSelectionStats, getDefaultChannelColors } from '../utilities/VivUtils.js';
import { ScatterplotLayer } from 'deck.gl';


class VivViewer {
  constructor(canvas,config,initialView){
    console.log('new VivViewer', config);
    this.canvas = canvas;


    
    this.height= this.canvas.height;
    this.width= this.canvas.width;
    this.config=config;
    loadOmeTiff(config.url,{pool:false}).then(loader=>{
      this._setUp(loader,initialView);
    });
  }

  setSize(x,y,conf){
    this.height=y;
    this.width=x;
    const v =this.getViewState(conf.x_scale,conf.y_scale,conf.offset);
    this.deck.setProps({
      height:y,
      width:x,
      viewState:v
    })
  }


  setPanZoom(offset,x_scale,y_scale){  
    const v= this.getViewState(x_scale,y_scale,offset);
    this.deck.setProps({
      viewState:v
    })
  }

  getViewState(x_scale,y_scale,offset){
    const hzoom = Math.log2(y_scale);
    const wzoom = Math.log2(x_scale);
    let xpos = ((1/x_scale)*(this.native_x))/2;
    xpos-= offset[0];
    let ypos = ((1/y_scale)*(this.native_y))/2;
    ypos+=this.native_y-offset[1];
    return {
      height:this.native_y,
      width:this.native_x,
      id:DETAIL_VIEW_ID,
      target:[xpos,ypos,0],
      zoom:[wzoom,hzoom]
    }
  }

  setChannel(channel){
    const channels=   this.layers[0].props;
    for (let i=0;i<channels.selections.length;i++){
      const sel= channels.selections[i];
      if (sel.c===channel.index){
        channels.colors[i]=hexToRGB(channel.color);
        channels.contrastLimits[i]=channel.contrastLimits;
        channels.channelsVisible[i]= channel.channelsVisible;
        break;
      }
    }
    this.layers=[this.layers[0]];
    this.deck.setProps({
      layers:[this.layers]
    })
  }

  removeChannel(channel){
    const chs=  this.layers[0].props;
    let i=0;
    for (let sel of chs.selections){
      if (sel.c===channel.index){
        break;
      }
      i++;
    }
    chs.colors.splice(i,1);
    chs.selections.splice(i,1);
    chs.contrastLimits.splice(i,1);
    chs.channelsVisible.splice(i,1);
    this.createLayers(chs);
    this.deck.setProps({
      layers:[this.layers]
    });

  }

  addChannel(channel){
    const chs=   this.layers[0].props;
    chs.channelsVisible.push(true);
    channel.color= channel.color || "#ff00ff";
    // pjt consider using helpers (copy from avivator utils).
    channel.contrastLimits = channel.contrastLimit || [20,100];
    channel.channelsVisible=true;
    chs.colors.push(hexToRGB(channel.color));
    chs.contrastLimits.push(channel.contrastLimits);
    chs.selections.push({z:0,t:0,c:channel.index});
    
    this.createLayers(chs);
    this.deck.setProps({
      layers:[this.layers],

    });

    channel.name=this.channels[channel.index].Name;
    return channel;

  }

  _setUpVolumeView(tiff) {
    const {SizeX, SizeY, SizeZ, Channels: channels} = tiff.metadata.Pixels;
    const target = [SizeX/2, SizeY/2, SizeZ/2];
    const id = '3d_' + DETAIL_VIEW_ID;
    const loader = tiff.data;
    const n = channels.length;
    const selections = channels.map((_, i) => {return {c: i, t: 0, z: 0}});
    const dtype = tiff.data[0].dtype;
    let { domains, contrastLimits } = getDefaultSelectionStats(n);
    getMultiSelectionStats(loader, selections).then((v) => {
      domains = v.domains;
      contrastLimits = v.contrastLimits;
      //TODO updateProps() equivalent
    });
    const colors = getDefaultChannelColors(n); //channels.map((_, i) => [i/n*255, (1-i/n)*255, 0]);
    const xSlice = [0, SizeX * 2];
    const ySlice = [0, SizeY * 2];
    const zSlice = [0, SizeZ * 2];
    const channelsVisible = channels.map(_ => true);
    const resolution = loader.length - 1;
    const extensions = [new ColorPalette3DExtensions.AdditiveBlendExtension()];
    const volumeView = new VolumeView({
      id,
      target,
      useFixedAxis: false,
      extensions,
      // extensions: [get3DExtension("", RENDERING_MODES.ADDITIVE)]
    });
    this.detailView = volumeView;
    const props = {
      id,
      loader,
      dtype,
      resolution,
      channelsVisible,
      contrastLimits,
      domains,
      selections,
      colors,
      xSlice, ySlice, zSlice
    };
    console.dir(props);
    const layers = volumeView.getLayers({
      props
    });
    this.volumeLayer = layers;
    this.volViewState = {
      zoom: 1, target
    };
    // const r = (v) => v * Math.random();
    // layers.push(new ScatterplotLayer({
    //   data: new Array(1000).fill().map(()=>{return {position: [r(SizeX), r(SizeY), r(SizeZ)]}}),
    //   radiusScale: 1,
    //   billboard: true,
    //   getFillColor: () => [100, 100, 100]
    // }))
  };

  _setUp(loader, iv){
    this.native_x= loader.metadata.Pixels.SizeX;
    this.native_y= loader.metadata.Pixels.SizeY;
    const {use3d} = this.config;
    
    this.extensions = [new ColorPaletteExtension()];
    this.channels = loader.metadata.Pixels.Channels;
    this.loader= loader.data;
    this.transparentColor=[255,255,255,0];
    const baseViewState = use3d ? undefined : this.getViewState(iv.x_scale,iv.y_scale,iv.offset);
    
    if (use3d) {
      this._setUpVolumeView(loader);
    } else {
      this.detailView = new DetailView({
        id: DETAIL_VIEW_ID,
        height:this.native_y,
        width:this.native_x
      });
    }
    const initialViewState = this.volViewState;
    const {image_properties} = this.config;
 
    const deckGLView =this.detailView.getDeckGlView();
    this.createLayers(image_properties);
    this.deck=new Deck({
          canvas:this.canvas,
          layers:[this.layers],
          views:[deckGLView],
          viewState:baseViewState,
          width:this.width,
          height:this.height,
          useDevicePixels:false,
          initialViewState,
          controller: use3d
    });
  }


  createLayers(info){
    //for now, pending refactor
    if (this.config.use3d) {
      this.layers = this.volumeLayer;
      return;
    }
    const viewStates=  {id: DETAIL_VIEW_ID}
    const layerConfig = {
      loader:this.loader,
      contrastLimits:info.contrastLimits.slice(0),
      colors:info.colors.slice(0),
      channelsVisible:info.channelsVisible.slice(0),
      selections:info.selections.slice(0),
      extensions:this.extensions,
      transparentColor:this.transparentColor
    };
    this.layers= this.detailView.getLayers({
      viewStates,
      props:layerConfig
    })

  }
}

export default VivViewer;