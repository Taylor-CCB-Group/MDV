const path = require('path');
webpack=require("webpack");

module.exports = {
        mode:"development",
		
        devtool: 'inline-source-map',
		devServer: {
			headers: {
				"Cross-Origin-Opener-Policy":"same-origin",
				"Cross-Origin-Embedder-Policy":"require-corp"		
			}
		  },
	
	   entry: './src/dev_index.js',

  	   output: {
    		   filename: 'dist/ciview2.js',

  	   },
  	   module:{ 
			
     		rules:[
			{
				test:/\.css$/,
				use:[
					"style-loader",
					"css-loader"
				]
				
			},
			{
				test: /\.(png|svg|jpg|gif|eot|ttf|woff|woff2)$/,
				type: 'asset/resource'
			
			}		
		]
  	}
};


