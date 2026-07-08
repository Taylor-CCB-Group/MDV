import reactHooks from "eslint-plugin-react-hooks";
import tseslint from "typescript-eslint";

const reactCompilerRules = reactHooks.configs["recommended-latest"].rules;

export default [
    {
        linterOptions: {
            reportUnusedDisableDirectives: "off",
        },
    },
    {
        ignores: [
            "node_modules/**",
            "dist/**",
            "vite-dist/**",
            "python/**",
            "src/react/components/avivatorish/**",
        ],
    },
    {
        files: ["src/**/*.{ts,tsx}"],
        languageOptions: {
            parser: tseslint.parser,
            parserOptions: {
                ecmaFeatures: {
                    jsx: true,
                },
                sourceType: "module",
            },
        },
        plugins: {
            "react-hooks": reactHooks,
        },
        rules: reactCompilerRules,
    },
];
