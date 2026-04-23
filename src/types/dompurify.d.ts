declare module "dompurify" {
    type PurifyProfiles = {
        svg?: boolean;
        svgFilters?: boolean;
        [key: string]: boolean | undefined;
    };

    type PurifyConfig = {
        USE_PROFILES?: PurifyProfiles;
        [key: string]: unknown;
    };

    type DOMPurify = {
        sanitize: (value: string, config?: PurifyConfig) => string;
    };

    const dompurify: DOMPurify;
    export default dompurify;
}
