import { action } from "mobx";
import { useChannelStats } from "../hooks";
import { VivMDVReact } from "./VivMDVReact";
import { observer } from "mobx-react-lite";
import { OmeTiffProvider, useChart, useOmeTiff } from "../context";
import { ChannelsState, VivProvider, useChannelsState, useChannelsStore, useImageSettingsStore } from "./avivatorish/state";


export default observer(function MainVivColorDialog() {
    return (
        <VivProvider>
            <Test />
        </VivProvider>
    )
})

export const Test = () => {
    const colors = useChannelsStore(({ colors }) => colors);
    const str = JSON.stringify(colors);
    return <div>{str}</div>
}