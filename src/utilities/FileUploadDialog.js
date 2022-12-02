import {BaseDialog} from "./Dialog.js";
import {createEl} from "../utilities/Elements.js";



 /**
  * A dialog which allows user to upload files 
  * @param {object} config - The settings for the dialog
  * @param {string} [config.title] - The dialog's title, if missing a generic title will be used 
  * @param {string} [config.filter] - The type of files allowed
  * @param {string} [config.text] Text describing which files to upload - a generic message will be
  * shown if absent
  * @param {function} [config.validate] - A function which accepts a list of files and validates
  * them - it should return an object with, files- a list of accepted files, msg - a string explaining
  * why files were invalid. If absent all files will be accepted
  * @param {function} [config.onstartupload] - A function which accepts the files to be uplaoded
  * and returns an object containing any data to send to the server
  * @param {function} config.url - The url for the upload
  */
class FileUploadDialog extends BaseDialog{
    constructor(config){
        const conf={
            footer:true,
            buttons:[{
                text:"Upload",
                method:"uploadfiles",
                id:"upload"
            }],
            title:config.title || "Upload Files"
        }
        const content={
           
            text:config.text || "Please press select to upload files"
        }
        super(conf,content);

        this.validate=config.validate;
        this.onupload= config.onupload;
        this.onfinished= config.onfinished;
        this.fileInput = createEl("input",{
            type:"file",
            multiple:true,
            filter:config.filter || "",
            styles:{},
            display:"none"
        },this.div);
        this.fileInput.addEventListener("change",e=>{
            this.filesChosen(e.target.files)
        });
        this.url = config.url;
        

    }

    init(content){
        this.files=[]
        createEl("div",{text:content.text,classes:["civ-dg-div"]},this.dialog);
        const b =createEl("div",{classes:["civ-dg-div","ciview-center"]},this.dialog);
        const s= createEl("button",{classes:["ciview-button-sm"],text:"Select Files"},b);
        createEl("i",{classes:["fas","fa-file-image"],styles:{marginLeft:"4px"}},s)
        b.addEventListener("click",()=>{
            this.fileInput.click()
        });
        this.fileDiv =createEl("div",{classes:["civ-dg-div"]},this.dialog);
        this.infoDiv= createEl("div",{classes:["civ-dg-div","ciview-center"]
        },this.dialog);
        this.progressDiv= createEl("progress",{classes:["civ-dg-div"],styles:{width:"100%"},value:0},this.dialog);
        this.disableButton("upload",true);
        
    }

    filesChosen(files){
        if (this.validate){
            const resp = this.validate(files);
            files = resp.files;
            this.infoDiv.textContent=resp.msg
        }
        for (let file of files){
            this.addFile(file);
        }

    }


    uploadfiles(){
        if (this.files.length>0){
            let data={}
            if (this.onupload){
                data = this.onupload(this.files);
            }
            let xhr = new XMLHttpRequest();
    		let fd = new FormData();
    		xhr.responseType="json";
    		xhr.open("POST", this.url, true);
            for (let index in  this.files){
				fd.append("file"+index, this.files[index]);
			}
    		fd.append("data",JSON.stringify(data));
    		xhr.onreadystatechange = ()=> {
        	    if (xhr.readyState == 4 && xhr.status == 200) {
                    let msg = null;
                    if (this.onfinished){
                        msg = this.onfinished(xhr.response);
                    }
                    this.infoDiv.innerHTML="";
                    msg = msg || "Files were successfully uploaded"
                    this.infoDiv.textContent=msg;
             }
            }
            xhr.upload.onprogress= e=>{
                this.progressDiv.value=e.loaded;
                this.progressDiv.max=e.total;
            }

            this.fileDiv.classList.add("ciview-disabled")

            this.disableButton("upload",true);
            this.infoDiv.innerHTML="";
            createEl("i",{
                classes:["fas","fa-spinner","fa-spin","ciview-center"],
                styles:{
                    fontSize:"30px"
                }
            },this.infoDiv)
           
    		xhr.send(fd);   

        }
        
    }

    removeFile(file){
        this.files=this.files.filter(x=>x!==file);
        if (this.files.length===0){
            this.disableButton("upload",true);
        }
    }

    addFile(file){
        if (this.files.find(x=>x.name===file.name)){
            return;
        }
        this.disableButton("upload",false);
        this.files.push(file);
        const h= createEl("div",{styles:{padding:"3px"}},this.fileDiv);
        createEl("span",{text:file.name},h);
        createEl("i",{
            classes:["fas","fa-trash","ciview-chart-icon"],
            styles:{float:"right"}
        },h)
        .addEventListener("click",()=>{
            h.remove();
            this.removeFile(file);
          
        })


    }


}

export default FileUploadDialog;