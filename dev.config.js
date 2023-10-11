const path = require('path');
webpack=require("webpack");
module.exports = (env) => {
	return {
        mode:"development",
		
        devtool: 'inline-source-map',
		devServer: {
			headers: {
				"Cross-Origin-Opener-Policy":"same-origin",
				"Cross-Origin-Embedder-Policy":"require-corp"
			},
			static:env.folder
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
  	}
}
};