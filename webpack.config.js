// PJT: for review: this seemed to be unused, was pointing at non-existent dev_index.js
// now making it equivalent to old dev.config.js,
// so devserver can be started without --config argument.

const path = require('path');
webpack=require("webpack");

module.exports = env =>{
	return{
        mode:"development",
		
        devtool: 'inline-source-map',
		devServer: {
			headers: {
				"Cross-Origin-Opener-Policy":"same-origin",
				"Cross-Origin-Embedder-Policy":"require-corp"		
			},
			static: env.folder || "./public/"
		  },
	
	   entry: './src/dev/devindex.js',

  	   output: {
    		   filename: 'dist/ciview2.js',

  	   },
		resolve: {
			extensions: ['.ts', '.js', '...']
		},
		module: {

			rules: [
				{
					test: /\.css$/,
					use: [
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