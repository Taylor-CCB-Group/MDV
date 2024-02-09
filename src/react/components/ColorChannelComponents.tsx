import { action } from "mobx";
import { observer } from "mobx-react-lite";
import { VivProvider, useChannelsStore, useChannelsStoreApi, useMetadata } from "./avivatorish/state";
import { Slider } from "@mui/material";


export default observer(function MainVivColorDialog() {
    return (
        <VivProvider>
            <Test />
        </VivProvider>
    )
})


const ChannelSliders = ({ index }: { index: number }) => {
    const limits = useChannelsStore(({ contrastLimits }) => contrastLimits);
    const {colors, selections} = useChannelsStore(({ colors, selections }) => ({colors, selections}));
    const metadata = useMetadata();
    if (!metadata) throw 'no metadata'; //TODO type metadata
    const name = metadata?.Pixels.Channels[selections[index].c].Name;
    const color = colors[index];
    const colorString = `rgb(${color[0]}, ${color[1]}, ${color[2]})`;
    const channelsStore = useChannelsStoreApi();
    // not sure I want to be using material.ui... consider adding a widget abstration layer.
    // not jumping right in with using it for all layout etc because I don't want to be tied to it.
    return (
        <div
        style = {{
            padding: '10px',
        }}
        >
            {name}
            <Slider
                style={{ color: colorString }}
                value={limits[index]}
                onChange={action((e, v) => {
                    limits[index] = v as [number, number];
                    const contrastLimits = [...limits];
                    channelsStore.setState({ contrastLimits });
                })}
            />
        </div>
    )
}

export const Test = () => {
    const colors = useChannelsStore(({ colors }) => colors);
    const colorsApi = useChannelsStoreApi();
    const str = JSON.stringify(colors);
    return <div style={{width: '100%'}}>{
        colors.map((c, i) => (
            <ChannelSliders key={i} index={i} />
        ))}
    </div>
}