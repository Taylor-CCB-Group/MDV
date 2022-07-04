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
			static:"./dist/hyperion_local/"
		  },
	
	   entry: './src/modules/hyperion_index.js',//'./_tests/tracks/browser.js',

  	   output: {
    		   filename: 'hyperion.js',
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
}
};