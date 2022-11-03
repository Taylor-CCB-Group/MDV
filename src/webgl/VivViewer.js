import {
  loadOmeTiff,
  DetailView,
  VolumeView,
  ColorPaletteExtension,
  ColorPalette3DExtensions,
  DETAIL_VIEW_ID,
} from '@hms-dbmi/viv';

import {hexToRGB, RGBToHex} from "../datastore/DataStore.js";
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
    this.initClip();
    loadOmeTiff(config.url).then(loader=>{
      this.tiff = loader;
      this._setUp(loader,initialView);
    });
  }

  setSize(x,y,conf){
    this.height=y;
    this.width=x;
    const v =this.getViewState(conf.x_scale,conf.y_scale,conf.offset);
    this.canvas.width = x;
    this.canvas.height = y;
    this.canvas.style.width = x;
    this.canvas.style.height = y;
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
    // when rendering 3d, we want viewState to be undefined so it can use initialViewState & internal camera control
    if (this.config.use3d) return undefined;
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

  /** equivalent to VivScatterPlot... */
  getAllChannels() {
    return this.channels;
  }
  getChannels() {
    const {props} = this.layers[0];
    const names = props.selections.map(x => this.channels[x.c].Name);
    const colors = props.colors.map(RGBToHex);
    return names.map((name, i) => {
      return {
        name,
        index: props.selections[i].c,
        color: colors[i],
        contrastLimits: props.contrastLimits[i].slice(0),
        channelsVisible: props.channelsVisible[i]
      }
    });
  }

  recenterCamera() {
    if (!this.config.use3d) return;
    console.log('recenter');
    const {SizeX, SizeY, SizeZ} = this.tiff.metadata.Pixels;
    const target = [SizeX / 2, SizeY / 2, SizeZ / 2];
    const initialViewState = {
      target,
      zoom: 1,
      rotationX: 0,
      rotationOrbit: 0 + Math.random() * 0.01
    }
    this.volViewState = initialViewState;
    this.deck.setProps({
      initialViewState
    });
  }

  _createLayers3D() {
    const tiff = this.tiff;
    //most of this can move into createLayers()
    const {SizeX, SizeY, SizeZ, Channels: channels} = tiff.metadata.Pixels;
    const target = [SizeX/2, SizeY/2, SizeZ/2];
    const id = '3d_' + DETAIL_VIEW_ID;
    const loader = tiff.data;
    const n = channels.length;
    const selections = channels.map((_, i) => {return {c: i, t: 0, z: 0}});
    const dtype = tiff.data[0].dtype;
    let { domains, contrastLimits } = getDefaultSelectionStats(n);
    ///wrong flow for now
    // getMultiSelectionStats(loader, selections).then((v) => {
    //   domains = v.domains;
    //   contrastLimits = v.contrastLimits;
    //   //TODO updateProps() equivalent
    // });
    const colors = getDefaultChannelColors(n); // this should change...
    const xSlice = this.getXSlice();
    const ySlice = this.getYSlice();
    const zSlice = this.getZSlice();
    const channelsVisible = channels.map(_ => true); // this should change...
    const resolution = loader.length - 1; // this should change...
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
    const volumeView = this.detailView;
    const layers = volumeView.getLayers({
      props
    });
    this.layers = layers;
    // const r = (v) => v * Math.random();
    // layers.push(new ScatterplotLayer({
    //   data: new Array(1000).fill().map(()=>{return {position: [r(SizeX), r(SizeY), r(SizeZ)]}}),
    //   radiusScale: 1,
    //   billboard: true,
    //   getFillColor: () => [100, 100, 100]
    // }))
  };

  initClip() {
    this.clipX = [0, 1];
    this.clipY = [0, 1];
    this.clipZ = [0, 1];
  }
  // so much boilerplate...
  setClipX(min, max) {
    this.clipX = [min, max];
    this._updateProps();
  }
  setClipY(min, max) {
    this.clipY = [min, max];
    this._updateProps();
  }
  setClipZ(min, max) {
    this.clipZ = [min, max];
    this._updateProps();
  }
  getXSlice() {
    const {SizeX} = this.tiff.metadata.Pixels;
    const [min, max] = this.clipX;
    const v = SizeX;
    return [min*v, max*v];
  }
  getYSlice() {
    const {SizeY} = this.tiff.metadata.Pixels;
    const [min, max] = this.clipY;
    const v = SizeY;
    return [min*v, max*v];
  }
  getZSlice() {
    const {SizeZ} = this.tiff.metadata.Pixels;
    const [min, max] = this.clipZ;
    const v = SizeZ;
    return [min*v, max*v];
  }
  _updateProps() {
    this.createLayers(); //as of this writing, this will lose changes made to channels etc.
    this.deck.setProps({layers: this.layers})
  }

  _setUp(loader, iv){
    this.native_x= loader.metadata.Pixels.SizeX;
    this.native_y= loader.metadata.Pixels.SizeY;
    const {use3d} = this.config;
    
    this.extensions = [new ColorPaletteExtension()];
    this.channels = loader.metadata.Pixels.Channels;
    this.loader= loader.data;
    this.transparentColor=[255,255,255,0];
    const baseViewState = this.getViewState(iv.x_scale,iv.y_scale,iv.offset);
    
    if (use3d) {
      // this._setUpVolumeView(loader);
      const { SizeX, SizeY, SizeZ } = loader.metadata.Pixels;
      const target = [SizeX / 2, SizeY / 2, SizeZ / 2];
      this.volViewState = {
        zoom: 1, target
      };
      this.detailView = new VolumeView({
        id: DETAIL_VIEW_ID,
        useFixedAxis: false,
        target,
        extensions: [new ColorPalette3DExtensions.AdditiveBlendExtension()],
      });
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
    this.createLayers(image_properties, loader);
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
    if (this.config.use3d) {
      this._createLayers3D();
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