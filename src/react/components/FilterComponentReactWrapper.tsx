import { g } from "@/lib/utils";
import { Autocomplete, TextField } from "@mui/material";
import { observer } from "mobx-react-lite";

export type DropdownType = {
    options: string[];
};

const FilterDropdown = observer(() => {
    const cm = window.mdv.chartManager;
    const { viewManager } = cm;
    const options = viewManager.all_views;
    return (
        <Autocomplete
        options={options}
        value={viewManager.current_view || null}
        onChange={(_event, newValue) => {
            if (newValue) {
              cm.changeView(newValue);
            }
          }}
        renderInput={(params) => <TextField {...params} label="Select View" />}
        sx={{display: "inline-flex", width: "20vw",
            margin: "0.2em",}}
        />
    );
});

const FilterComponentReactWrapper = () => {
    return (
            <FilterDropdown />
        // <div style={{maxWidth: '10em', display: 'inline-flex'}}>
        // {/* <DropdownAutocompleteComponent props={g({
        //     label: '',
        //     current_value: "default",
        //     type: 'dropdown',
        //     values: values,
        // })} /> */}
        // </div>
    );
};

export default FilterComponentReactWrapper;