import type { ColumnSelectionProps, CTypes } from "@/lib/columnTypeHelpers";
import { observer } from "mobx-react-lite";
import { useRowsAsColumnsLinks } from "../chartLinkHooks";
import { useCallback, useEffect, useMemo, useState } from "react";
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
    const { link, linkedDs } = props.link;
    const subgroupKeys = useMemo(
        () => Object.keys(link.subgroups),
        [link.subgroups],
    );
    const firstSubgroupKey = subgroupKeys[0] ?? "";

    const defaultMaxItems = props.multiple ? 10 : 1;
    const [queryObj] = useState(() => {
        if (props.current_value) {
            const v = isArray(props.current_value)
                ? props.current_value[0]
                : props.current_value;
            if (v instanceof RowsAsColsQuery && v.link === link) return v;
        }
        return new RowsAsColsQuery(link, linkedDs.name, defaultMaxItems);
    });

    /** Subgroup tab selection; applied onto `queryObj` only when the user clicks Activate (like column/link mode). */
    const [pendingSubgroupKey, setPendingSubgroupKey] = useState(
        () => queryObj.effectiveSubgroupKey,
    );

    // Keep the tab in sync when the chart's active query for this link changes from outside.
    useEffect(() => {
        const ov = props.current_value;
        const v = isArray(ov) ? ov[0] : ov;
        if (v instanceof RowsAsColsQuery && v.link === link) {
            setPendingSubgroupKey(v.effectiveSubgroupKey);
        }
    }, [props.current_value, link]);

    const isActiveResolved = useMemo(() => {
        const ov = props.current_value;
        const v = isArray(ov) ? ov[0] : ov;
        return (
            typeof v !== "string" &&
            v instanceof RowsAsColsQuery &&
            v.link === link &&
            v.effectiveSubgroupKey === pendingSubgroupKey &&
            v.maxItems === queryObj.maxItems
        );
    }, [
        props.current_value,
        link,
        pendingSubgroupKey,
        queryObj.maxItems,
    ]);

    const setMaxItems = useCallback(
        action((v: number) => {
            queryObj.maxItems = v;
        }),
        [queryObj],
    );

    const activateLink = useCallback(() => {
        action(() => {
            queryObj.subgroupName =
                pendingSubgroupKey === firstSubgroupKey ? "" : pendingSubgroupKey;
        })();
        //@ts-expect-error setSelectedColumn generic not behaving nicely
        props.setSelectedColumn(props.multiple ? [queryObj] : queryObj);
    }, [
        props,
        queryObj,
        pendingSubgroupKey,
        firstSubgroupKey,
    ]);

    const subgroupLabels = useMemo(
        () =>
            subgroupKeys.map(
                (k) => link.subgroups[k]?.label ?? k,
            ),
        [subgroupKeys, link.subgroups],
    );

    const pendingLabel =
        link.subgroups[pendingSubgroupKey]?.label ?? pendingSubgroupKey;

    return {
        isActive: isActiveResolved,
        activateLink,
        maxItems: queryObj.maxItems,
        setMaxItems,
        link,
        queryObj,
        subgroupKeys,
        subgroupLabels,
        pendingSubgroupKey,
        setPendingSubgroupKey,
        pendingLabel,
    };
}

const ActiveLinkComponent = observer(
    <T extends CTypes, M extends boolean>(props: LocalProps<T, M>) => {
        const {
            isActive,
            activateLink,
            maxItems,
            setMaxItems,
            subgroupKeys,
            subgroupLabels,
            pendingSubgroupKey,
            setPendingSubgroupKey,
            pendingLabel,
            queryObj,
        } = useActiveLink(props);

        return (
            <div className="w-full flex flex-col gap-1 text-xs font-light p-1">
                {subgroupKeys.length > 1 ? (
                    <TabHeader
                        activeTab={pendingSubgroupKey}
                        setActiveTab={setPendingSubgroupKey}
                        tabs={subgroupKeys}
                        tabLabels={subgroupLabels}
                        indicatorTab={queryObj.effectiveSubgroupKey}
                    />
                ) : null}
                <div className="flex justify-end items-center gap-1">
                    <Button onClick={activateLink}>
                        {isActive ? "Using" : "Activate"} &apos;
                        {pendingLabel}&apos; link
                    </Button>
                    {props.multiple && (
                        <TextField
                            size="small"
                            className="max-w-20 float-right"
                            type="number"
                            label="Max"
                            variant="standard"
                            value={maxItems}
                            onChange={(e) =>
                                setMaxItems(Number(e.target.value))
                            }
                        />
                    )}
                </div>
            </div>
        );
    },
);

const ActiveLinkMultiComponent = observer(
    <T extends CTypes, M extends boolean>(
        props: ColumnSelectionProps<T, M>,
    ) => {
        const links = useRowsAsColumnsLinks();
        const committedLinkName = useMemo(() => {
            const v = props.current_value;
            const cur = isArray(v) ? v[0] : v;
            if (cur instanceof RowsAsColsQuery) {
                return cur.link.name;
            }
            return undefined;
        }, [props.current_value]);

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
                <TabHeader
                    activeTab={activeLinkName}
                    setActiveTab={handleLinkChange}
                    tabs={links.map((l) => l.link.name)}
                    indicatorTab={committedLinkName}
                />

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
