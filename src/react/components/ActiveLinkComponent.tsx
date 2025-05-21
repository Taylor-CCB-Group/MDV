import type { ColumnSelectionProps, CTypes } from "@/lib/columnTypeHelpers";
import { observer } from "mobx-react-lite";
import { useRowsAsColumnsLinks } from "../chartLinkHooks";
import { useCallback, useMemo, useState } from "react";
import { isArray } from "@/lib/utils";
import { RowsAsColsQuery } from "@/links/link_utils";
import { Button, TextField } from "@mui/material";
import { action } from "mobx";
import TabHeader from "./TabHeader";

type LocalProps<T extends CTypes, M extends boolean> = ColumnSelectionProps<
    T,
    M
> & {
    link: ReturnType<typeof useRowsAsColumnsLinks>[0];
};

function useActiveLink<T extends CTypes, M extends boolean>(
    props: LocalProps<T, M>,
) {
    const { link, linkedDs } = props.link; //pending multiple link support
    const sg = Object.values(link.subgroups)[0]; //pending subgroup support
    // const highlightedForeignRows = useHighlightedForeignRows();
    const isActive = useMemo(() => {
        const ov = props.current_value;
        const v = isArray(ov) ? ov[0] : ov;
        // still needs more understanding of sg etc etc
        return (
            typeof v !== "string" &&
            v instanceof RowsAsColsQuery &&
            v.link === link
        ); // && v.subgroup === sg;
    }, [props.current_value, link]);
    // if we already have a link, we should initialize the maxItems to its current value
    const defaultMaxItems = props.multiple ? 10 : 1;
    // will change when we can actually select different link/sg... but we could just keep one of these around
    // apply it to the app state / mutate whenever... this is very un-functional/react... which is a shame (really).
    const [queryObj] = useState(() => {
        if (props.current_value) {
            const v = isArray(props.current_value)
                ? props.current_value[0]
                : props.current_value;
            if (v instanceof RowsAsColsQuery && v.link === link) return v;
        }
        return new RowsAsColsQuery(link, linkedDs.name, defaultMaxItems);
    });
    const setMaxItems = useCallback(
        action((v: number) => {
            queryObj.maxItems = v;
        }),
        [],
    );
    const activateLink = useCallback(() => {
        //@ts-expect-error setSelectedColumn generic not behaving nicely
        props.setSelectedColumn(props.multiple ? [queryObj] : queryObj);
    }, [props, queryObj]);

    // todo
    // const freezeLink = useCallback(() => {
    // })
    return {
        isActive,
        activateLink,
        maxItems: queryObj.maxItems,
        setMaxItems,
        link,
        sg,
        // highlightedForeignRows,
    };
}

const ActiveLinkComponent = observer(
    <T extends CTypes, M extends boolean>(props: LocalProps<T, M>) => {
        const {
            isActive,
            activateLink,
            maxItems,
            setMaxItems,
            sg,
            // highlightedForeignRows,
        } = useActiveLink(props);
        // maybe something like rowsText can be shown in a little scrollable box...
        // const rowsText = highlightedForeignRows.slice(0, maxItems).map((r) => r.name).join(", ");

        return (
            <div className="w-full flex justify-end text-xs font-light p-1">
                <Button onClick={activateLink}>
                    {isActive ? "Using" : "Activate"} '{sg.label}' Link
                </Button>
                {props.multiple && (
                    <TextField
                        size="small"
                        className="max-w-20 float-right"
                        type="number"
                        label="Max"
                        variant="standard"
                        value={maxItems}
                        onChange={(e) => setMaxItems(Number(e.target.value))}
                    />
                )}
                {/* {isActive && rowsText} */}
            </div>
        );
    },
);

const ActiveLinkMultiComponent = observer(
    <T extends CTypes, M extends boolean>(
        props: ColumnSelectionProps<T, M>,
    ) => {
        const links = useRowsAsColumnsLinks();
        // Infer the initial active tab based on `current_state`
        const [activeLinkName, setActiveLinkName] = useState(() => {
            const v = props.current_value;
            const currentFieldSpec = isArray(v) ? v[0] : v;
            if (currentFieldSpec instanceof RowsAsColsQuery) {
                return currentFieldSpec.link.name;
            }
            return links[0].link.name;
        });

        const handleLinkChange = useCallback((name: string) => {
            setActiveLinkName(name);
        }, []);

        // Don't show tab UI if we only have one link
        if (links.length === 1) {
            return (
                <ActiveLinkComponent
                    {...props}
                    link={links[0]}
                />
            );
        }

        return (
            <div>
                {/* Tab-based interface for selecting the active link */}
                <TabHeader
                    activeTab={activeLinkName}
                    setActiveTab={handleLinkChange}
                    tabs={links.map((link) => link.link.name)}
                />

                {/* Render the active link's UI */}
                {links.map((link) =>
                    link.link.name === activeLinkName ? (
                        <ActiveLinkComponent
                            key={link.link.name}
                            {...props}
                            link={link}
                        />
                    ) : null
                )}
            </div>
        );
    },
);

export default ActiveLinkMultiComponent;
