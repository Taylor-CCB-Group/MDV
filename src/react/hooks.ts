import { useEffect, useMemo, useState } from "react";
import BaseChart from "../charts/BaseChart";

export function useChartSize(chart: BaseChart) {
    const div = chart.contentDiv;
    const [width, setWidth] = useState(div.clientWidth);
    const [height, setHeight] = useState(div.clientHeight);
    useEffect(() => {
        const resize = () => {
            setWidth(div.clientWidth);
            setHeight(div.clientHeight);
        };
        const observer = new ResizeObserver(resize);
        observer.observe(div);
        return () => observer.unobserve(div);
    }, [div]);
    return [width, height];
}

