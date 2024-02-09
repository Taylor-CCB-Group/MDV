import { action } from "mobx";
import { observer } from "mobx-react-lite";
import { VivProvider, useChannelsStore, useChannelsStoreApi } from "./avivatorish/state";
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
    const channelsStore = useChannelsStoreApi();
    // not sure I want to be using material.ui... consider adding a widget abstration layer.
    return (
        <div>
            <Slider
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