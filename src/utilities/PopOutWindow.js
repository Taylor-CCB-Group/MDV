/**
* A popout window
@deprecated
*/
class PopOutWindow {
    /**
    * @param {function} onload - the function called when the window has loaded passing the
    * window's document and and object conatining the width and height of the window.
    * @param {function} onclose - the function called when the window closes passing the
    * window's document.
    * @param {object} config - setup information for the window.
    * @param {int} [config.url] - By default a url is created using createObjectURL, however
    * some browsers prevent scripts from running in the window as it is considered cross domain.
    * You can supply a url (pointing to an empty/small html page) from the local doamin to overcome this
    * e.g '/pages/popout.html'
    * @param {int} [config.width=800] - the width of the window.
    * @param {int} [config.height=400] - the height of the window.
    * @param {int} [config.screenX=200] - the x position of the window.
    * @param {int} [config.screenY=200] - the y position of the window.
    * @param {function} [config.onresize] - a function that is called when the window is resized.
    * The function has the following parameters:- the document element and an object containing
    * the current width and height of the window.
    * @param {function} [config.stayInFront=false] - if true then the window will always stay in front of 
    * the parent, however as it temporarily loses focus, its tends to flicker.  
    * @param {string} [config.title] - the title of the window.
    * @param {string[]} [config.stylesheets] - a list of style sheet locations to add to the header.
    * (any styles/stylesheets in the parent window are are automatically copied)
    */

    constructor(onload,onclose,config={}) {
        let self = this;
        const className = document.documentElement.className;
        // let winHtml = `<!DOCTYPE html><html class="${className}"></html>`; //not used
        let winUrl = URL.createObjectURL(new Blob([],{
            type: "text/html"
        }));
        winUrl = config.url || winUrl;
   
        const c= config;

        this.window = window.open(winUrl, "win" + new Date().getTime(), 
            `width=${c.width || 800},
            height=${c.height || 400},
            screenX=${c.screenX||200},
            screenY=${c.screenY||200}
            ;modal=yes`);

        this.window.addEventListener("load",function() {
            let doc = self.window.document;
            doc.documentElement.className = className; //something unsets this again... not sure what.
            // doc.body.className = className; //also gets unset.
            this.setTimeout(()=> {
                doc.documentElement.className = className;
            }, 1000);
            doc.write(`<head>`);
            //copy header from parent window
            let h = document.head.innerHTML;
            doc.write(h);
            //add any additional style sheets
            if (c.stylesheets){
                for (let sheet of c.stylesheets) {
                    doc.write(`<link rel="stylesheet" type="text/css" href="${sheet}">`);
                }
            }
            doc.write(`</head>`);   
            doc.write(`<body><body>`);
            if (c.title){
                doc.title=c.title;
            }
      
            if (onload){
                onload(doc,self.getDimensions());
            } 
            self.onclose=onclose;    
            self.window.addEventListener("beforeunload",function(e) {
                    self.closing(doc);
            });
            
            if (c.stayInFront){
                self.focusListener= ()=>{
                    self.window.focus();
                };
                document.body.addEventListener("mousedown",self.focusListener)
              
            }
            const resized = c.onresize;
           
            if (resized){
                self.window.addEventListener("resize",function(e) {
                    //resizing is usually expensive so throttle it 
                    //to stop the callback being called frequently
                    clearTimeout(self.resizeTO);
                    self.resizeTO=setTimeout(()=>{
                       
                        resized(doc,self.getDimensions());
                    },200);
                  
                });     
            }         
        });
    }

    getDimensions(){
        const w = this.window;
        return {width:w.innerWidth-10,height:w.innerHeight-10};
    }

    closing(doc){
        if (this.onclose){
            this.onclose(doc);
        }
        if (this.focusListener){
            document.body.removeEventListener("mousedown",this.focusListener);
        }
    }

    /**
    * closes the window
    */
    close(){
        this.window.close();
    }
}
export {PopOutWindow}