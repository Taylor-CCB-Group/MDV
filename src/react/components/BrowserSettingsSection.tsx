import { type BrowserSetting, browserSettings } from "@/lib/browserSettings";
import { FormControl, FormHelperText, InputLabel, MenuItem, Select, Typography } from "@mui/material";
import { observer } from "mobx-react-lite";

const ChoiceSettingControl = observer(function ChoiceSettingControl({ setting }: { setting: BrowserSetting }) {
    const labelId = `browser-setting-${setting.key}`;
    return (
        <FormControl size="small" fullWidth sx={{ mt: 2 }}>
            <InputLabel id={labelId}>{setting.label}</InputLabel>
            <Select
                labelId={labelId}
                label={setting.label}
                value={setting.get()}
                onChange={(event) => setting.set(event.target.value)}
            >
                {setting.choices.map((choice) => (
                    <MenuItem key={choice.value} value={choice.value}>
                        {choice.label}
                    </MenuItem>
                ))}
            </Select>
            {setting.description && <FormHelperText>{setting.description}</FormHelperText>}
        </FormControl>
    );
});

/** Every setting in `browserSettings`; changes apply at once and are kept in this browser. */
export default function BrowserSettingsSection() {
    return (
        <>
            <Typography>
                These apply to every project opened in this browser, and are not saved with projects or views.
            </Typography>
            {browserSettings.map((setting) => (
                <ChoiceSettingControl key={setting.key} setting={setting} />
            ))}
        </>
    );
}
