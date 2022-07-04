const path = require('path');
webpack=require("webpack");

module.exports = {	  
	   entry: './src/modules/basic_index.js',	 
  	   output: {
    		   path: path.resolve(__dirname,"dist/basic"),
    		   filename: 'ciview2.js',
			   assetModuleFilename: 'images/[name][ext]',
			  
  	   },	  
  	   module:{ 			
     		rules:[
			{
				test:/\.css$/,
				use:[
					"style-loader",
					"css-loader"
				]
			} ,
			{
				test: /\.(png|svg|jpg|gif|eot|ttf|woff|woff2)$/,
				type: 'asset/resource',
				generator:{
					publicPath:"./"
				}		
			  }		
		]		
  	}
	 
};
