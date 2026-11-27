import { LegendControlButton } from "@/react/components/legend/LegendFrame";

/** The buttons on the strip of legends above a chart's plot. */
export default function LegendStripControls({ onFloat }: { onFloat: () => void }) {
    return <LegendControlButton icon="fas fa-clone" label="Float legends over the plot" onClick={onFloat} />;
}
