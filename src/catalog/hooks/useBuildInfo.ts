export type BuildInfo = {
    commitDate: string;
    branchName: string;
    commitHash: string;
    lastCommitMessage: string;
    buildDate: string;
    dirty: boolean;
}
const useBuildInfo = () => {
    /**
     * Returns information about the current build, particularly git revision.
     * 
     * Note, hardcoded for vite, this will need to be updated if we switch to a different bundler.
     */
    function getBuildInfo(): BuildInfo {
        const commitDate = import.meta.env.VITE_GIT_COMMIT_DATE || "";
        const branchName = import.meta.env.VITE_GIT_BRANCH_NAME || "";
        const commitHash = import.meta.env.VITE_GIT_COMMIT_HASH || "";
        const lastCommitMessage = import.meta.env.VITE_GIT_LAST_COMMIT_MESSAGE || "";
        const buildDate = import.meta.env.VITE_BUILD_DATE || "";
        const dirty = import.meta.env.VITE_GIT_DIRTY === "dirty";
        return {
            commitDate,
            branchName,
            commitHash,
            lastCommitMessage,
            buildDate,
            dirty,
        };
    }

    return { buildInfo: getBuildInfo() };
};

export default useBuildInfo;


