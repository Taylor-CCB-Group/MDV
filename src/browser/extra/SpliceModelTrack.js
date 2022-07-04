
import { MLVBedTrack,MLVTrack,Graphics } from "../tracks";

class SpliceModelTrack extends MLVBedTrack{
	constructor(config){
        config.featureHeight= config.featureHeight || 18;
		super(config);
		this._setFeatureSource();
		this.filter_function=null;
		this.color_function=null;
		
	}

   

    renderFeature(feature, coord,ctx,info){
      
        if (info.bpPerPixel>100) {
            if (this.color_function){
                ctx.fillStyle=this.color_function(feature);
            }
            ctx.fillRect(coord.px, coord.py, coord.pw, coord.h);
            //Graphics.fillRect(ctx, p0, y_pos1, p1 - p0, 10, {fillStyle: c});
        }
        else {
            const data = this.dataStore.getRowAsObject(feature.data[0],["snp_alt","donor_pos_gain","donor_pos_loss","acceptor_pos_gain","acceptor_pos_loss"])
            let snp_pos = Math.floor((parseInt(feature.data[1])+1 - info.bpStart) / info.bpPerPixel);
            const w = Math.round(1 / info.bpPerPixel);
            const w2= w<3?3:w;
            
            let don_gain = Math.round((data.donor_pos_gain -1- info.bpStart) / info.bpPerPixel);
            ctx.fillStyle = 'green';
            ctx.fillRect(don_gain, coord.py+6, w2, coord.h-6);
            don_gain = Math.round((data.donor_pos_loss - info.bpStart) / info.bpPerPixel);
            ctx.fillStyle = 'red';
            ctx.fillRect(don_gain, coord.py+6, w2, coord.h-6);
            don_gain = Math.round((data.acceptor_pos_gain - info.bpStart) / info.bpPerPixel);
            ctx.fillStyle = '#7FFF00';
            ctx.fillRect(don_gain, coord.py+6, w2, coord.h-6);
            don_gain = Math.round((data.acceptor_pos_loss - info.bpStart) / info.bpPerPixel);
            ctx.fillStyle = 'orange';
            ctx.fillRect(don_gain, coord.py+6, w2, coord.h-6);
            snp_pos-=w/2;

            Graphics.strokeText(ctx, data.snp_alt, snp_pos, coord.py, {
                strokeStyle: Graphics.nucleotideColors[data.snp_alt],
                font: 'normal 12px Arial',
                textAlign: 'center',
                textBaseline:'top'
            });
        }  
    }
}


MLVTrack.track_types["splice_model_track"]={
	"class":SpliceModelTrack,
	name:"Splce Model"
}