import JsonView from "react18-json-view";
import useDiff from "./utils/useDiff";

export type StateDiffProps = {
    input1: object;
    input2: object;
};

const StateDiffComponent = ({ input1, input2 }: StateDiffProps) => {

    const diff = useDiff(input1, input2);

    return (
        <>
            <JsonView
                style={{
                    backgroundColor: "transparent",
                    fontSize: "0.875rem",
                    fontFamily: "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace",
                }}
                src={diff}
                collapseStringsAfterLength={10}
                collapsed={1}
            />
        </>
    );
};

export default StateDiffComponent;
