import { action } from "mobx";
import { useChannelStats } from "../hooks";
import { VivMDVReact } from "./VivMDVReact";
import { observer } from "mobx-react-lite";
import { OmeTiffProvider, useChart, useOmeTiff } from "../context";
import { ChannelsState, useChannelsState } from "../viv_state";
import { Button, InputLabel, MenuItem, Select, Slider } from "@mui/material";


export default observer(function MainVivColorDialog() {
    return (
        <OmeTiffProvider>
            <div>
                <ChannelContrastEditor index={0} />
                <ChannelSelect />
            </div>
        </OmeTiffProvider>
    )
})


const ChannelSelect = observer(() => {
    const { config } = useChart() as VivMDVReact;
    if (config.type === 'VivMdvRegionReact') {
        return <div>Color channel selection not available for region views.</div>
    }
    else if (config.type !== 'VivMdvReact') throw new Error('unexpected config type');
    const ome = useOmeTiff();
    if (!ome) return <div>loading...</div>;
    const channelOptions = ome.metadata.Pixels.Channels.map((c, i) => (
        <MenuItem key={c.ID} value={i}>{c.Name}</MenuItem>
    ));
    return (
    <div>
        <InputLabel id="channel-select-label">Channel</InputLabel>
        <Select value={config.channel} labelId="channel-select-label" label="Channel" onChange={
            action(e => config.channel = Number.parseInt(e.target.value))}>
            {channelOptions}
        </Select>
        <Button onClick={action(() => {
            // add channel...
        })}>Add channel</Button>
    </div>
    )
});

const ChannelContrastEditor = observer(function ChannelContrastEditor({index}: {index: number}) {
    const channelsState = useChannelsState();
    const ome = useOmeTiff();
    const stats = useChannelStats(ome, 0);
    //nb, this is unused copilot code... and the 'index' is not used sensibly...
    const value = channelsState.contrastLimits[index] || stats?.contrastLimits;
    return (
        value && <Slider
            getAriaLabel={() => `Channel ${index} contrast`}
            value={value}
            min={0}
            max={65536}
            onChange={action((e, v) => channelsState.contrastLimits[index] = v as [number, number])}
        />
    )
});
