import { useEffect, useMemo, useState } from "react";
import BaseChart from "../charts/BaseChart";
import { createRoot } from "react-dom/client";

let ID = 0;
function ReactTest({dataStore}) {
  const [id] = useState(ID++);
  const [filterSize, setFilterSize] = useState(dataStore.filterSize);
  useEffect(() => {
    const listener = () => setFilterSize(dataStore.filterSize);
    dataStore.addListener(id, listener);
    return () => dataStore.removeListener(id, listener);
  }, [dataStore]);
  return (
    <div>
      <h1>React Test</h1>
      <h2>filterSize: {filterSize}</h2>
    </div>
  );
}

class ReactChart extends BaseChart {
  constructor(dataStore, div, config) {
    super(dataStore, div, config);
    createRoot(this.contentDiv).render(<ReactTest dataStore={dataStore} />);
  }
}

BaseChart.types["ReactChart"] = {
  "class": ReactChart,
  name: "React HelloWorld Chart",
}

// we rely on the side-effect of this import to register the chart type
export default 42;