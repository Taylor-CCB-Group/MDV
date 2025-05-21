import type { AxiosResponse } from "axios";

export type ParseErrorProps = {
    response: Response | AxiosResponse;
    fallbackText: string;
}

export const parseErrorResponse = async ({response, fallbackText}: ParseErrorProps) => {
    let errorData;
    if (response instanceof Response) {
        errorData = await response.json().catch(() => ({
            error: fallbackText,
        }));
    } else {
        errorData = response.data;
    }

    if (response.status === 403) {
        return new Error(errorData?.error || "Forbidden: You are not allowed to perform this action.");
    } else if (response.status === 404) {
        return new Error(errorData?.error || "Resource Not Found.");
    } else if (response.status === 500) {
        return new Error(errorData?.error || "Internal Server Error.");
    } else {
        return new Error(errorData?.error || fallbackText);
    }
};