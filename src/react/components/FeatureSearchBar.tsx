import ClearIcon from "@mui/icons-material/Clear";
import SearchIcon from "@mui/icons-material/Search";
import { Chip, IconButton, InputAdornment, TextField, Tooltip, Typography } from "@mui/material";
import { action } from "mobx";
import { observer } from "mobx-react-lite";
import type { ClipboardEvent, KeyboardEvent } from "react";
import type { DataType, LoadedDataColumn } from "../../charts/charts";
import { useChart } from "../context";
import { isSearchActive, mergeTerms, parseFeatureList } from "./featureSearch";
import { clearFeatureSearch, useFeatureSearchFilter, useMissingTerms } from "./featureSearchHooks";
import type { GeneNetworkConfig } from "./GeneNetworkChart";

const MAX_CHIPS = 30;

/**
 * Search box for the Gene Info chart. Typing filters the view by substring; pasting a list,
 * or pressing Enter, adds exact identifiers that the view is then limited to.
 */
export const FeatureSearchBar = observer(function FeatureSearchBar({
    column,
    shownCount,
}: {
    column: LoadedDataColumn<DataType>;
    shownCount: number;
}) {
    const { config } = useChart<GeneNetworkConfig>();
    const search = config.search;
    const terms = search.terms ?? [];
    useFeatureSearchFilter(column);
    const missing = useMissingTerms(column, terms);

    const addTerms = action((added: string[]) => {
        search.terms = mergeTerms(terms, added);
    });
    const setQuery = action((query: string) => {
        search.query = query;
    });
    const removeTerm = action((term: string) => {
        search.terms = terms.filter((t) => t !== term);
    });

    const handlePaste = (event: ClipboardEvent<HTMLInputElement>) => {
        const items = parseFeatureList(event.clipboardData.getData("text"));
        // a single word is just typing; a list becomes terms
        if (items.length < 2) return;
        event.preventDefault();
        addTerms(items);
        // a leftover query would also have to match, hiding the pasted terms
        setQuery("");
    };

    const handleKeyDown = (event: KeyboardEvent<HTMLInputElement>) => {
        if (event.key === "Enter") {
            const items = parseFeatureList(search.query ?? "");
            if (items.length === 0) return;
            event.preventDefault();
            addTerms(items);
            setQuery("");
        } else if (event.key === "Backspace" && !search.query && terms.length > 0) {
            removeTerm(terms[terms.length - 1]);
        }
    };

    return (
        <div className="flex flex-col gap-1 px-2 pt-2 pb-1">
            <TextField
                size="small"
                fullWidth
                value={search.query ?? ""}
                onChange={(e) => setQuery(e.target.value)}
                onPaste={handlePaste}
                onKeyDown={handleKeyDown}
                placeholder="Search, or paste a list"
                slotProps={{
                    htmlInput: { "aria-label": "Search identifiers" },
                    input: {
                        startAdornment: (
                            <InputAdornment position="start">
                                <SearchIcon fontSize="small" />
                            </InputAdornment>
                        ),
                        endAdornment: isSearchActive(search) ? (
                            <InputAdornment position="end">
                                <Tooltip title="Clear search">
                                    <IconButton size="small" onClick={() => clearFeatureSearch(config)}>
                                        <ClearIcon fontSize="small" />
                                    </IconButton>
                                </Tooltip>
                            </InputAdornment>
                        ) : undefined,
                    },
                }}
            />
            {terms.length > 0 && (
                <div className="flex flex-wrap gap-1 max-h-20 overflow-y-auto">
                    {terms.slice(0, MAX_CHIPS).map((term) => (
                        <Chip
                            key={term}
                            size="small"
                            label={term}
                            color={missing.includes(term) ? "warning" : "default"}
                            variant={missing.includes(term) ? "outlined" : "filled"}
                            onDelete={() => removeTerm(term)}
                        />
                    ))}
                    {terms.length > MAX_CHIPS && (
                        <Chip size="small" variant="outlined" label={`+${terms.length - MAX_CHIPS} more`} />
                    )}
                </div>
            )}
            {isSearchActive(search) && (
                <Typography variant="caption" color="text.secondary">
                    {shownCount} shown
                    {missing.length > 0 && (
                        <Tooltip title={missing.join(", ")}>
                            <Typography component="span" variant="caption" color="warning.main">
                                {` · ${missing.length} not found`}
                            </Typography>
                        </Tooltip>
                    )}
                </Typography>
            )}
        </div>
    );
});
