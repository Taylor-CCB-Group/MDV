import { styled } from '@mui/material/styles';
import { TextField, type TextFieldProps } from "@mui/material";

const StyledAdornment = styled('div')(({ theme }) => ({
    display: 'flex',
    alignItems: 'center',
    // border: `1px solid ${theme.palette.divider}`,
    // boxSizing: 'border-box',
    '.MuiAutocomplete-endAdornment &': {
        marginRight: theme.spacing(1),  // Adjust spacing as needed
    }
}));

/** Modified version of TextField that allows a `customEndAdornment`
 * along with the standard endAdornment passed in `InputProps`.
 */
export const TextFieldExtended = (props: TextFieldProps & { customEndAdornment?: React.JSX.Element; customStartAdornment?: React.JSX.Element}) => {
    const { InputProps, customEndAdornment, ...rest } = props;
    const inputProps: typeof InputProps = {
        ...InputProps,
        endAdornment: (
            <StyledAdornment>{customEndAdornment} {InputProps?.endAdornment}</StyledAdornment>
        ),
        startAdornment: (
            <StyledAdornment>{props.customStartAdornment} {InputProps?.startAdornment}</StyledAdornment>
        )
    };
    return <TextField {...rest} InputProps={inputProps} />;
};
