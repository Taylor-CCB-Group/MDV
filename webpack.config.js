// PJT: for review: this seemed to be unused, was pointing at non-existent dev_index.js
// now making it equivalent to old dev.config.js,
// so devserver can be started without --config argument.

const path = require('path');
webpack=require("webpack");

module.exports = env =>{
	const conf  = {
		mode:"development",	
		entry: "./src/modules/desktop_index.js",

		output: {
			filename: 'static/js/mdv.js',
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
	if (env.dev){
		const folder = env.dev.startsWith("folder")?env.dev.replace("folder:",""):null;
		conf.devtool= 'inline-source-map';
		conf.devServer= {
			headers: {
				"Cross-Origin-Opener-Policy":"same-origin",
				"Cross-Origin-Embedder-Policy":"require-corp"		
			},
			static:folder?folder:"./public/"
		  }
		  if (!folder){
			conf.devServer.proxy= [
				{
					context: ['/images', '/tracks', '/get_configs', '/get_view',
					 '/get_binary_data', '/get_row_data', '/get_data', '/save_state',],
					target: env.dev,
				}
			]
		}
	}
	else{
		conf.mode=env.mode || "production";
		if (env.build==="desktop"){
			conf.output= {
				path: path.resolve(__dirname,"python/mdvtools/static/js"),
				filename: 'mdv.js',
				assetModuleFilename: '../img/[name][ext]',	   
			}
		}
		else if (env.build==="production"){
			conf.entry= env.nofont?'./src/modules/basic_index_nf.js':'./src/modules/basic_index.js';
			conf.output= {
				path: env.outputpath || path.resolve(__dirname, "dist/basic"),
				filename: 'mdv.js',
				assetModuleFilename: 'images/[name][ext]',
			};
		//where the js looks for assets - default is images folder in same location as page
			conf.module.rules[1].generator={
				publicPath:env.assetpath || "./"
			}
		}
	}
	return conf;
};