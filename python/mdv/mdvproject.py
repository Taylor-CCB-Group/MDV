import os
import h5py
import numpy
import json
import gzip
from os.path import join,split
from  shutil import copytree,ignore_patterns
from .server import create_app

class MDVProject:
    def __init__(self,dir):
        self.dir=dir
        self.h5file = join(dir,"datafile.h5")
        self.datasourcesfile= join(dir,"datasources.json")
        self.statefile= join(dir,"state.json")
        self.viewsfile= join(dir,"views.json")
        self.imagefolder = join(dir,"images")
        self.trackfolder = join(dir,"tracks")
        self._datasources=None
        #if os.path.exists(dir):
        #    os.mkdir(dir)

    @property
    def datasources(self):
        if not self._datasources:
            if not os.path.exists(self.datasourcesfile):
                self._datasources=[]
            else:
                self._datasources=get_json(self.datasourcesfile)
        return self._datasources

    @datasources.setter
    def datasources(self,value):
        self._datasources=value
        save_json(self.datasourcesfile,value)



    def get_links(self,datasource,filter=None):
        ds =  [x for x in self.datasources if x["name"]==datasource][0]
        links=[]
        lnks= ds.get("links")
        if lnks:
            for lnkto in lnks:
                lnk = lnks[lnkto]
                if (filter==None or lnk.get(filter)):
                    links.append({
                        "datasource":lnkto,
                        "link":lnk
                    })
        return links
                    

    def serve(self,**kwargs):
        create_app(self,**kwargs)

    
    def get_configs(self):
        config ={
            "datasources":get_json(self.datasourcesfile),
            "state":get_json(self.statefile)
        }
        #legacy 
        hyperion_conf= join(self.dir,"hyperion_config.json")
        if os.path.exists(hyperion_conf):
            config["hyperion_config"]= get_json(hyperion_conf)
        #end
        return config

    
    def get_view(self,view):
            views =get_json(self.viewsfile)
            return views[view]
    

    def convert_to_static_page(self,outdir):   
        copytree(self.dir,outdir,ignore=ignore_patterns("*.h5"))
        self.convert_data_to_binary(outdir)
        template = join(split(os.path.abspath(__file__))[0],"templates","page.html")
        page = open(template).read()
        page =page.replace('_mdvInit()','_mdvInit(true)')
        ifile = join(outdir,"index.html")
        o = open(ifile,"w")
        o.write(page)
        o.close()




    def convert_data_to_binary(self,outdir=None):
        if not outdir:
            outdir=self.dir
        h5 =  h5py.File(self.h5file)
        dss = self.datasources
        for ds in dss:
            n = ds["name"]
            gr = h5[n]
            dfile = join(outdir,"{}.b".format(n))
            o = open(dfile,"wb")
            index={}
            current_pos=0
            for c in ds["columns"]:     
                dt = gr.get(c["field"])
                if not dt:
                    continue
                arr = numpy.array(dt)
                comp = gzip.compress(arr.tobytes())
                o.write(comp)
                new_pos = current_pos +len(comp)
                index[c["field"]]=[current_pos,new_pos-1]
                current_pos = new_pos

            #add rows to columns gene score / tensors etc
            lnks = self.get_links(n,"rows_as_columns")
            for ln in lnks:
                rc= ln["link"]["rows_as_columns"]
                for sg in  rc["subgroups"]:
                    info = rc["subgroups"][sg]
                    sgrp = gr[info["name"]]
                    sparse = info.get("type")=="sparse"
                    #get number of rows in linked datasource
                    plen = [x["size"] for x in dss if x["name"]==ln["datasource"]][0]
                    for i in range (0,plen):
                        comp=   gzip.compress(get_subgroup_bytes(sgrp,i,sparse))
                        o.write(comp)
                        new_pos = current_pos +len(comp)
                        index[f'{sg}{i}']=[current_pos,new_pos-1]
                        current_pos = new_pos
                  
            o.close()    
            ifile = dfile[:dfile.rindex(".")]+".json"
            i = open (ifile,"w")
            i.write(json.dumps(index))
            i.close()

    def get_byte_data(self,columns,group):
        h5 = h5py.File(self.h5file,"r")
        byte_list=[]  
        for column in columns:
            sg = column.get("subgroup")
            if sg:
                sgindex= int(column["sgindex"])
                byte_list.append(get_subgroup_bytes(h5[group][sg],sgindex,column.get("sgtype")=="sparse"))
            else:
                data = h5[group][column["field"]]      
                byte_list.append(numpy.array(data).tobytes())         
            
        h5.close()
        return b''.join(byte_list)


def get_json(file):
    return json.loads(open(file).read())

def save_json(file,data):
    o = open(file,"w")
    o.write(json.dumps(data,indent=2))
    o.close()

def get_subgroup_bytes(grp,index,sparse=False):
  
    if sparse:
        offset = grp["p"][index:index+2]
        _len = offset[1]-offset[0]
        _indexes = numpy.array(grp["i"][offset[0]:offset[1]])
        _values=  numpy.array(grp["x"][offset[0]:offset[1]],numpy.float32)
        return numpy.array([_len],numpy.uint32).tobytes()  \
                        + numpy.array(_indexes).tobytes() \
                        +  numpy.array(_values).tobytes()
    else:
        _len =grp["length"][0]
        offset= index*_len
        return numpy.array(grp["x"][offset:offset+_len],numpy.float32).tobytes()
