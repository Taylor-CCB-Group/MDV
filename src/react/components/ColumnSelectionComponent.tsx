import { useMemo, useState } from "react";
import { observer } from "mobx-react-lite";
import TextField from "@mui/material/TextField";
import Autocomplete from "@mui/material/Autocomplete";
import { useDataStore } from "../context.js";

type ColumnSelectionProps = {
    setSelectedColumn: (column: string) => void;
    placeholder?: string;
    exclude?: string[];
};
/**
 * A component for selecting a column from the data store.
 * Must be in a context where `useDataStore` is available
 * (e.g. if we're in a more global dialog etc rather than a chart context,
 * this would be ambiguous).
 */
export default observer((props: ColumnSelectionProps) => {
    const { setSelectedColumn, placeholder } = props;
    const dataStore = useDataStore();
    // todo column groups
    const columns = useMemo(
        () => dataStore.columns
            .filter((c) => !props.exclude?.includes(c.name))
            // .filter((c) => c.datatype === "multitext") //todo datatype prop etc
            .map((c) => c.name),
        [dataStore, props.exclude],
    );
    return (
        <>
            <Autocomplete
                options={columns}
                onChange={(_, value) => {
                    setSelectedColumn(value);
                }}
                renderInput={(params) => {
                    const { key, ...p } = params as typeof params & {
                        key: string;
                    };
                    return (
                        <TextField
                            key={key}
                            {...p}
                            // label="Annotation column"
                            placeholder={placeholder}
                        />
                    );
                }}
                renderOption={(props, text) => {
                    const { key, ...p } = props as typeof props & {
                        key: string;
                    };
                    return (
                        <li key={key} {...p}>
                            {text}
                        </li>
                    );
                }}
            />
        </>
    );
});
