import BaseChart from "./BaseChart";
import { createEl } from "../utilities/Elements.js";

const applySelectableTextStyle = (element) => {
    element.style.userSelect = "text";
    element.style["-webkit-user-select"] = "text";
    element.style["-moz-user-select"] = "text";
    element.style["-ms-user-select"] = "text";
};

class TextBoxChart extends BaseChart {
    constructor(dataStore, div, config) {
        super(dataStore, div, config);
        const c = this.config;
        c.text = c.text || "";
        this.para = createEl(
            "div",
            {
                styles: {
                    padding: "5px",
                    overflowY: "auto",
                    maxHeight: "100%",
                },
            },
            this.contentDiv,
        );
        //conditional import of function which uses the packages marked and sanitize
        //decreases size of entry module
        import("../utilities/MarkdownText").then(({ default: renderText }) => {
            this.render = renderText;
            this.textContainer = createEl(
                "div",
                {},
                this.para,
            )
            if (!this.textContainer) return;
            this.textContainer.innerHTML = this.render(c.text);
            applySelectableTextStyle(this.textContainer.firstChild);

        });
    }

    getSettings() {
        const settings = super.getSettings();
        const c = this.config;
        settings.push({
            label: "Markdown Text",
            type: "textbox",
            current_value: c.text,
            func: (x) => {
                c.text = x;
                this.para.innerHTML = this.render(x);
            },
        });
        return settings;
    }
}

export default TextBoxChart;

BaseChart.types["text_box_chart"] = {
    class: TextBoxChart,
    name: "Text Box",
    params: [],
    extra_controls: () => [
        {
            type: "textbox",
            name: "text",
            label: "Markdown Text:",
        },
    ],
    init: (config, dataSource, extraControls) => {
        config.text = extraControls.text;
    },
};
