import { FormControl, FormHelperText, InputLabel, MenuItem, Select, Typography } from "@mui/material";
import { observer } from "mobx-react-lite";
import { type BrowserSetting, browserSettings } from "@/lib/browserSettings";
import { matchesQuery } from "./debugDialogSearch";

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

/** Whether a setting matches a search query, by label, description, key or choice labels. */
export function browserSettingMatchesQuery(setting: BrowserSetting, query: string) {
    return matchesQuery(query, setting.label, setting.description, setting.key, ...setting.choices.map((c) => c.label));
}

/** The given settings (default: all of `browserSettings`); changes apply at once and are kept in this browser. */
export default function BrowserSettingsSection({
    settings = browserSettings,
}: {
    settings?: readonly BrowserSetting[];
}) {
    return (
        <>
            <Typography>
                These apply to every project opened in this browser, and are not saved with projects or views.
            </Typography>
            {settings.map((setting) => (
                <ChoiceSettingControl key={setting.key} setting={setting} />
            ))}
        </>
    );
}
