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
                "Cross-Origin-Embedder-Policy": "require-corp"
            },
            static: "./examples/"
        },

        entry: './src/modules/obvios_index.ts',

        output: {
            filename: 'dist/ciview2.js',

        },
        plugins: [
            // new BundleAnalyzerPlugin()
        ],
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