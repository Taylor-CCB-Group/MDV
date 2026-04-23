import { observer } from "mobx-react-lite";
import { useChart } from "@/react/context";
import type { TextBoxChartConfig } from "../TextBoxChartReactWrapper";
import TextBoxMarkdownRenderer from "./TextBoxMarkdownRenderer";

const TextBoxChartComponent = observer(() => {
    const chart = useChart<TextBoxChartConfig>();

    return (
        <div className="mdv-textbox-root" data-textbox-chart="true">
            <div className="mdv-textbox-scroll">
                <TextBoxMarkdownRenderer markdown={chart.config.text ?? ""} />
            </div>
        </div>
    );
});

export default TextBoxChartComponent;
