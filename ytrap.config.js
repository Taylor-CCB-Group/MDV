const path = require('path');
webpack = require("webpack");
const BundleAnalyzerPlugin = require('webpack-bundle-analyzer').BundleAnalyzerPlugin;


module.exports = env => {
    return {
        mode: "development",

        devtool: 'inline-source-map',
        devServer: {
            headers: {
                "Cross-Origin-Opener-Policy": "same-origin",
                "Cross-Origin-Embedder-Policy": "require-corp",
                "Access-Control-Allow-Origin": "*",
                "Access-Control-Allow-Methods": "GET",
                "Access-Control-Allow-Headers": "X-Requested-With, content-type, Authorization",
            },
            static: ["./examples/", "../ytrap2"],
        },

        entry: './src/modules/ytrap_static_index.ts',

        output: {
            filename: 'dist/ciview2.js',

        },
        plugins: [
            // new BundleAnalyzerPlugin()
        ],
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