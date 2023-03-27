const path = require('path');
webpack=require("webpack");

module.exports = (env)=> {
	return {	  
	   	entry: './src/modules/desktop_index.js',
  	   	output: {
    		   path: path.resolve(__dirname,"python/mdv/static/mdv/js"),
    		   filename: 'mdv.js',
			   assetModuleFilename: '../img/[name][ext]',
			  
  	   	},
		resolve: {
			extensions: ['.ts', '.js', '...']
		},
		mode:env.mode,  
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
						publicPath:"/static/mdv/img/"
					}		
				},	
				{
					test: /\.tsx?$/,
					use: 'ts-loader',
					exclude: '/node_modules'
				}	
			]		
  		}
	}
}