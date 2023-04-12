import { useEffect, useMemo, useState } from "react";
import BaseChart from "../charts/BaseChart";
import { createRoot } from "react-dom/client";

let ID = 0;
function ReactTest({dataStore, parent}) {
  const [id] = useState(ID++);
  const [filterSize, setFilterSize] = useState(dataStore.filterSize);
  const [text, setText] = useState(parent.config.text);
  useEffect(() => {
    const listener = () => setFilterSize(dataStore.filterSize);
    dataStore.addListener(id, listener);
    parent.addListener("text", (type, data) => setText(data));
    return () => {
      dataStore.removeListener(id, listener);
      parent.removeListener("text");
    }
  }, [dataStore]);
  return (
    <div>
      <h1>React Test</h1>
      <h2>filterSize: {filterSize}</h2>
      <p>{text}</p>
    </div>
  );
}

class ReactChart extends BaseChart {
  constructor(dataStore, div, config) {
    super(dataStore, div, config);
    createRoot(this.contentDiv).render(<ReactTest dataStore={dataStore} parent={this} />);
  }

  getSettings() {
    const c = this.config;
    const settings = super.getSettings();
    return settings.concat([
      {
        label:"Text",
        type:"textbox",
        current_value: c.text,
        func: x => {
          c.text = x;
          this._callListeners("text", x);
        }
      }
    ]);
  }
}

BaseChart.types["ReactChart"] = {
  "class": ReactChart,
  name: "React HelloWorld Chart",
  init: (config, dataSource, extraControls) => {
    config.text = "Hello World";    
  }
}

// we rely on the side-effect of this import to register the chart type
export default 42;