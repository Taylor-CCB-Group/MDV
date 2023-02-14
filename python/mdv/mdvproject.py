import os
import h5py
import numpy
import json
from mdv.server import create_app

class MDVProject:
    def __init__(self,dir):
        self.dir=dir
        self.h5file = os.path.join(dir,"datafile.h5")
        self.datsourcesfile= os.path.join(dir,"datasources.json")
        self.statefile= os.path.join(dir,"state.json")
        self.viewsfile= os.path.join(dir,"views.json")
        self.imagefolder = os.path.join(dir,"images")
        self.trackfolder = os.path.join(dir,"tracks")
        #if os.path.exists(dir):
        #    os.mkdir(dir)


    def serve(self,**kwargs):
        create_app(self,**kwargs)

    
    def get_configs(self):
        config ={
            "datasources":get_json(self.datsourcesfile),
            "state":get_json(self.statefile)
        }
        #legacy 
        hyperion_conf= os.path.join(self.dir,"hyperion_config.json")
        if os.path.exists(hyperion_conf):
            config["hyperion_config"]= get_json(hyperion_conf)
        #end
        return config

    
    def get_view(self,view):
            views =get_json(self.viewsfile)
            return views[view]
    

    def get_byte_data(self,columns,group):
        h5 = h5py.File(self.h5file,"r")
        byte_list=[]  
        for column in columns:
            sg = column.get("subgroup")
            if sg:
                sgindex= int(column["sgindex"])
                if column.get("sgtype")=="sparse":
                    offset = h5[group][sg]["p"][sgindex:sgindex+2]
                    _len = offset[1]-offset[0]
                    _indexes = numpy.array(h5[group][sg]["i"][offset[0]:offset[1]])
                    _values=  numpy.array(h5[group][sg]["x"][offset[0]:offset[1]],numpy.float32)
                    byte_list.append(numpy.array([_len],numpy.uint32).tobytes()  \
                        + numpy.array(_indexes).tobytes() \
                        +  numpy.array(_values).tobytes())
                else:
                    _len =h5[group][sg]["length"][0]
                    offset= sgindex*_len
                    byte_list.append(numpy.array(h5[group][sg]["x"][offset:offset+_len],numpy.float32).tobytes())
            else:
                data = h5[group][column["field"]]      
                byte_list.append(numpy.array(data).tobytes())         
            
        h5.close()
        return b''.join(byte_list)


def get_json(file):
    return json.loads(open(file).read())