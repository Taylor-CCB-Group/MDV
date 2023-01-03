const path = require('path');
webpack = require("webpack");
module.exports = (env) => {
	return {
		mode: "production",

		entry: './src/modules/hyperion_index.js',

		output: {
			path: path.resolve(__dirname, "dist/hyperion_local"),
			filename: 'hyperion.js',
			assetModuleFilename: 'hyperion_assets/[name][ext]',

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