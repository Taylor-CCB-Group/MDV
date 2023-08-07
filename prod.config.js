const path = require('path');
webpack = require("webpack");
module.exports = (env) => {
    return {
        entry: './src/modules/basic_index.js',
        output: {
            path: env.outputpath || path.resolve(__dirname, "dist/basic"),
            filename: 'mdv.js',
            assetModuleFilename: 'images/[name][ext]',
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
                    type: 'asset/resource',
                    generator: {
                        publicPath:env.assetpath || "./"
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
};
