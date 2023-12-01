const path = require('path');
webpack=require("webpack");
module.exports = (env) => {
	const config=  {
        mode:"development",
		
        devtool: 'inline-source-map',
		devServer: {
			headers: {
				"Cross-Origin-Opener-Policy":"same-origin",
				"Cross-Origin-Embedder-Policy":"require-corp"
			}
			//static:env.folder
		  },
	
	   entry: "./src/modules/desktop_index.js",

  	   output: {
    		   filename: 'static/js/mdv.js',
  	   },
		 resolve: {
            extensions: ['.ts', '.js', '...']
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
				
				},
				{
					test: /\.tsx?$/,
					use: 'ts-loader',
					exclude: '/node_modules'
				}
			]
		},
	}
	if (env.folder){
		config.devServer.static=env.folder;
	}
	if (env.server){
		config.devServer.proxy= [
			{
				context: ['/images', '/tracks', '/get_configs', '/get_view',
				 '/get_binary_data', '/get_row_data', '/get_data', '/save_state',],
				target: env.server,
			}
		]
	}
	return config
};