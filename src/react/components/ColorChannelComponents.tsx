import { useId, useState } from "react";
import { VivProvider, useChannelsStore, useChannelsStoreApi, useLoader, useMetadata } from "./avivatorish/state";
// some of this is quite bloated - could use dynamic imports for some of the more complex components
import { Checkbox, FormControl, InputLabel, MenuItem, Select, Slider } from "@mui/material";
import { PopoverPicker } from "./ColorPicker";

export default function MainVivColorDialog() {
    return (
        <VivProvider>
            <Test />
        </VivProvider>
    )
}


const ChannelChooserMUI = ({ index }: { index: number }) => {
    const channels = useMetadata().Pixels.Channels.map(c => c.Name) as string[];
    const selections = useChannelsStore(({selections}) => selections);
    const channelsStore = useChannelsStoreApi();
    const id = useId();
    const name = channels[selections[index].c]

    return (
        <>
        <FormControl fullWidth variant="standard">
            <InputLabel id={id}>{name}</InputLabel>
            <Select
            labelId={id}
            label="Channel"
            onChange={e => {
                const newSelections = [...selections];
                if (typeof e.target.value !== 'number') {
                    return;
                }
                newSelections[index].c = e.target.value;
                channelsStore.setState({selections: newSelections});
            }}
            >
                {channels.map((c, i) => (<MenuItem key={`${i}_${c}`} value={i}>{c}</MenuItem>))}
            </Select>
        </FormControl>
        </>
    )
}

const ChannelChooser = ({index}: {index: number}) => {
    const channels = useMetadata().Pixels.Channels.map(c => c.Name) as string[];
    const selections = useChannelsStore(({selections}) => selections);
    const channelsStore = useChannelsStoreApi();

    return (
        <>
        <select
        value={selections[index].c}
        style={{width: '100%', padding: '0.2em'}} onChange={e => {
            const newSelections = [...selections];
            try {
                newSelections[index].c = Number.parseInt(e.target.value)
                channelsStore.setState({selections: newSelections});
            } catch {}
        }}>
            {channels.map((c, i) => (<option key={`${i}_${c}`} value={i}>{c}</option>))}
        </select>
        </>
    )
}

const ChannelController = ({ index }: { index: number }) => {
    const limits = useChannelsStore(({ contrastLimits }) => contrastLimits);
    const {colors, domains, channelsVisible, removeChannel} = useChannelsStore(({ colors, domains, channelsVisible, removeChannel }) => (
        { colors, domains, channelsVisible, removeChannel }
    )); //channelOptions is a string[] of channel names...
    const metadata = useMetadata();
    const channelsStore = useChannelsStoreApi();

    if (!metadata) throw 'no metadata'; //TODO type metadata
    const channelVisible = channelsVisible[index];
    const color = colors[index];
    const colorString = `rgb(${color[0]}, ${color[1]}, ${color[2]})`;
    // not sure I want to be using material.ui... consider adding a widget abstration layer.
    // not jumping right in with using it for all layout etc because I don't want to be tied to it.
    // consider adding tailwind soon...
    return (
        <>
        <div
        style = {{
            // padding: '10px',
            display: 'grid',
            gridTemplateColumns: '0.4fr 0.1fr 0.1fr 1fr 0.1fr',
            justifyItems: 'flex-start',
            alignItems: 'center',
        }}
        >
            <ChannelChooser index={index} />
            <Checkbox
                checked={channelVisible}
                onChange={() => {
                    channelsVisible[index] = !channelVisible;
                    const visible = [...channelsVisible];
                    channelsStore.setState({ channelsVisible: visible });
                }}
            />
            <PopoverPicker color={color} onChange={c => {
                colors[index] = c;
                const newColors = [...colors];
                channelsStore.setState({ colors: newColors });
            }}
            />
            <Slider
                size="small"
                style={{ color: colorString, marginLeft: '10px' }}
                value={limits[index]}
                min={domains[index][0]}
                max={domains[index][1]}
                onChange={(_, v) => {
                    limits[index] = v as [number, number];
                    const contrastLimits = [...limits];
                    channelsStore.setState({ contrastLimits });
                }}
            />
            <button
            style={{marginLeft: '12px'}}
            onClick={() => {
                removeChannel(index);
            }}>
                x
            </button>
        </div>
        </>
    )
}

const AddChannel = () => {
    // const loader = useLoader();
    // const { labels } = loader[0];
    // const channelsStore = useChannelsStoreApi();
    const addChannel = useChannelsStore(state => state.addChannel);
    return <button onClick={() => {
        addChannel({
            selections: {c: 0, z: 0, t: 0},
            ids: String(Math.random()),
            channelsVisible: true,
            colors: [255, 255, 255],
        }); //todo isLoading, auto-contrast, palette, etc.
    }}>
        Add channel
    </button>;
}

export const Test = () => {
    const colors = useChannelsStore(({ colors }) => colors);
    return <div style={{width: '100%'}}>{
        colors.map((_, i) => (
            <ChannelController key={i} index={i} />
        ))}
        <AddChannel />
    </div>
}
