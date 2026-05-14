import type { ColumnSelectionProps, CTypes } from "@/lib/columnTypeHelpers";
import { observer } from "mobx-react-lite";
import type { DropdownMappedValues, GuiSpec } from "@/charts/charts";
import { useCallback, useEffect, useMemo, useState } from "react";
import { action, makeAutoObservable, observe, runInAction } from "mobx";
import { DropdownAutocompleteComponent } from "./SettingsDialogComponent";
import { g, isArray } from "@/lib/utils";
import { getFieldName, subgroupKeyFromFieldName, RowsAsColsQuery, type RowsAsColslink } from "@/links/link_utils";
import { useRowsAsColumnsLinks } from "../chartLinkHooks";
import TabHeader from "./TabHeader";

type LocalProps<T extends CTypes, M extends boolean> = ColumnSelectionProps<
    T,
    M
> & {
    link: ReturnType<typeof useRowsAsColumnsLinks>[0];
};

function pickInitialSubgroupKey(
    subgroups: RowsAsColslink["subgroups"],
    current_value: unknown,
): string {
    const keys = Object.keys(subgroups);
    const first = keys[0] ?? "";
    const v = isArray(current_value) ? current_value[0] : current_value;
    if (typeof v !== "string") {
        return first;
    }
    const prefix = subgroupKeyFromFieldName(v);
    if (prefix && subgroups[prefix]) {
        return prefix;
    }
    return first;
}

/**
 * A hook for managing the state of manually chosen columns from a rows-as-columns link.
 *
 * Should encapsulate
 * - making sure that an appropriate type of value is selected
 * - providing a way to select the value that has a simple interface, managing the complexity of incoming generic types
 *
 * Subgroup tabs only change which subgroup's **options** are listed in the dropdown; they do not change the chart
 * until the user picks a field (same idea as column / link / active-link modes). The saved `FieldName` strings
 * carry their own subgroup prefix; we validate against that prefix, not the tab alone.
 *
 * ! there is a bug - when values from a different link are selected, the spec doesn't update
 * we want to be able to set more different combinations of values - which will mean more changes to state management
 * for now, the UI should better reflect what the user is allowed to select.
 */
function useLinkSpec<T extends CTypes, M extends boolean>(props: LocalProps<T, M>) {
    const { linkedDs, link } = props.link;
    const { name_column, name, subgroups } = link;
    const targetColumn = linkedDs.dataStore.columnIndex[name_column];
    if (!targetColumn) {
        throw new Error(`Could not find name_column '${name_column}' in linked dataset`);
    }
    if (linkedDs.size > 2 ** 16) {
        // If we have more than 2^16 rows, I think that means there must be some duplicates...
        // so this is potentially somewhat undefined behaviour.
        console.warn("Linked dataset has more than 2^16 rows, this may lead to unexpected behaviour");
    }

    const subgroupKeys = useMemo(() => Object.keys(subgroups), [subgroups]);
    const subgroupLabels = useMemo(
        () => subgroupKeys.map((k) => subgroups[k]?.label ?? k),
        [subgroupKeys, subgroups],
    );

    /** Which subgroup's rows populate the dropdown (browse only until a field is chosen). */
    const [browseSubgroupKey, setBrowseSubgroupKey] = useState(() =>
        pickInitialSubgroupKey(subgroups, props.current_value),
    );

    /**
     * Subgroup prefix embedded in the chart's current `FieldName` value(s), if any.
     * We do not auto-sync `browseSubgroupKey` when this changes: that would snap the tab back right after the user
     * picks a feature on another layer (bound updates from the new field). Initial browse still uses
     * `pickInitialSubgroupKey`; use "Show …" below to align the tab with the saved selection.
     */
    const boundSubgroupKey = useMemo(() => {
        const v = props.current_value;
        const first = isArray(v) ? v[0] : v;
        if (typeof first !== "string" || !first.includes("|")) {
            return null;
        }
        const p = subgroupKeyFromFieldName(first);
        return p && subgroups[p] ? p : null;
    }, [props.current_value, subgroups]);

    const sg = browseSubgroupKey;

    const values = useMemo(() => {
        if (!link.valueToRowIndex) {
            // because in certain instances the link isn't ready yet...
            // we can handle that here, but ideally we should do it earlier in the process
            console.warn("Link not ready yet - waiting to get values for UI... should be handled differently.");
            link.initPromise.then(() => {
                console.log(link.valueToRowIndex);
            });
            return [[], "label", "fieldName"] as DropdownMappedValues<"label", "fieldName">;
        }
        const labelValueObjects = [...link.valueToRowIndex.entries()].map(([value, rowIndex]) => ({
            label: value,
            fieldName: getFieldName(sg, value, rowIndex),
        }));
        return [[...labelValueObjects], "label", "fieldName"] satisfies DropdownMappedValues<
            "label",
            "fieldName"
        >;
    }, [sg, link.valueToRowIndex, link.initPromise]);

    const { setSelectedColumn } = props;
    const isMultiType = props.multiple;

    const setValue = useCallback(
        (v: string | string[]) => {
            //@ts-expect-error setSelectedColumn generic not behaving nicely
            setSelectedColumn(v);
        },
        [setSelectedColumn],
    );

    const getSafeInternalValue = useCallback(
        (v: typeof props.current_value) => {
            const defaultVal = props.multiple ? [] : "";
            if (!v) return defaultVal;

            const rowValueOk = (field: string, prefix: string): boolean => {
                if (!field.startsWith(`${prefix}|`)) {
                    return false;
                }
                const rowVal = field.split("|")[1].split(`(${prefix})`)[0].trimEnd();
                return Boolean(link.valueToRowIndex?.has(rowVal));
            };

            if (props.multiple) {
                if (!isArray(v)) {
                    return defaultVal;
                }
                if (v.length === 0) {
                    return defaultVal;
                }
                const [firstValue] = v;
                if (typeof firstValue !== "string") {
                    return defaultVal;
                }
                if (!firstValue.includes("|")) {
                    return defaultVal;
                }
                const prefix = subgroupKeyFromFieldName(firstValue);
                if (!prefix || !link.subgroups[prefix]) {
                    return defaultVal;
                }
                for (const item of v) {
                    if (typeof item !== "string" || !rowValueOk(item, prefix)) {
                        return defaultVal;
                    }
                }
                return v;
            }
            if (isArray(v)) {
                return defaultVal;
            }
            if (typeof v !== "string") {
                return defaultVal;
            }
            if (!v.includes("|")) {
                return defaultVal;
            }
            const prefix = subgroupKeyFromFieldName(v);
            if (!prefix || !link.subgroups[prefix]) {
                return defaultVal;
            }
            if (!rowValueOk(v, prefix)) {
                return defaultVal;
            }
            return v;
        },
        [link.valueToRowIndex, link.subgroups, props.multiple],
    );

    const [initialValue] = useState(() => getSafeInternalValue(props.current_value));

    /// we don't want a new spec every time the value changes... we want to give it an observable value
    // (I think)- definition of this spec may want to be in the hook, with an observable value that isn't
    // the real props.current_value???
    const [spec] = useState(() =>
        makeAutoObservable(
            g<"multidropdown" | "dropdown">({
                type: isMultiType ? "multidropdown" : "dropdown",
                label: `specific '${name}'`, //don't really want to have this label here, would take a bigger change to remove
                values,
                current_value: initialValue,
                func: action((_v) => {
                    // something about having setValue as a dependency makes it end up re-running this...
                    // & we really want the spec to be stable.
                    //! for some reason if we call this in here we get glitchy behaviour, not really sure
                    // setValue(v);
                }),
            }),
        ),
    );

    useEffect(() => {
        runInAction(() => {
            spec.values = values;
        });
    }, [spec, values]);

    useEffect(() => {
        // do you swim with the mobx stream? or against it?
        // took me a while to realise that this is the direction we should be going in...
        return observe(spec, "current_value", (change) => {
            if (change.type === "update") {
                setValue(change.newValue);
            }
        });
    }, [spec, setValue]);

    const handleSubgroupTab = useCallback((key: string) => {
        setBrowseSubgroupKey(key);
    }, []);

    return {
        spec,
        showSubgroupTabs: subgroupKeys.length > 1,
        subgroupKeys,
        subgroupLabels,
        browseSubgroupKey,
        handleSubgroupTab,
        boundSubgroupKey,
    };
}

/**
 * This presents the user with a method for selecting column(s) corresponding to a linked dataset.
 *
 * There should be a re-usable component for selecting values from a category column, and this should be used here,
 * referring to the `link.name_column` column of the linked dataset... at the moment, we deal with this in a more ad-hoc way.
 *
 * Values that this component operates on - at the time of writing - will be in a `FieldName` (`string`) format, with pipes separating
 * `subgroup|value (subgroup)|index`, where `subgroup` is the name of the user-selected subgroup, `value (subgroup)` is a human-readable
 * representation of the value, and `index` is the index of the value in the column.
 */
const LinkToColumnComponent = observer(<T extends CTypes, M extends boolean>(props: LocalProps<T, M>) => {
    const {
        spec,
        showSubgroupTabs,
        subgroupKeys,
        subgroupLabels,
        browseSubgroupKey,
        handleSubgroupTab,
        boundSubgroupKey,
    } = useLinkSpec(props);
    return (
        <div className="flex flex-col p-1 gap-1">
            {showSubgroupTabs ? (
                <TabHeader
                    activeTab={browseSubgroupKey}
                    setActiveTab={handleSubgroupTab}
                    tabs={subgroupKeys}
                    tabLabels={subgroupLabels}
                    indicatorTab={boundSubgroupKey ?? undefined}
                />
            ) : null}
            <DropdownAutocompleteComponent
                props={spec as GuiSpec<"dropdown"> | GuiSpec<"multidropdown">}
            />
        </div>
    );
});

/**
 * Avoid dealing with internal state of the link until it's ready.
 */
const LinkToColumnComponentInit = observer(<T extends CTypes, M extends boolean>(props: LocalProps<T, M>) => {
    const { link } = props.link;
    //wanted to initialise this with ~init.linkPromise.resolved but apparently that's not a thing
    //so we have an initial re-render even if the link is already resolved.
    const [ready, setReady] = useState(false);
    useEffect(() => {
        link.initPromise.then(() => {
            setReady(true);
        });
    }, [link.initPromise]);
    if (!ready) {
        return null;
    }
    return <LinkToColumnComponent {...props} />;
});

const LinkMulti = observer(<T extends CTypes, M extends boolean>(props: ColumnSelectionProps<T, M>) => {
    const links = useRowsAsColumnsLinks();
    const committedLinkIndicator = useMemo(() => {
        const v = props.current_value;
        const cur = isArray(v) ? v[0] : v;
        if (cur instanceof RowsAsColsQuery) {
            return cur.link.name;
        }
        if (typeof cur === "string") {
            const prefix = subgroupKeyFromFieldName(cur);
            if (!prefix) {
                return undefined;
            }
            return links.find((entry) => Boolean(entry.link.subgroups[prefix]))?.link.name;
        }
        return undefined;
    }, [props.current_value, links]);

    const [activeLinkName, setActiveLinkName] = useState(() => {
        // Determine the initial active link based on `current_state` and `FieldName`
        // this logic might be adjusted if we allow combined selections with multiple links
        // also we might want to persist as user changes tab outside of this component
        const v = props.current_value;
        const currentFieldName = isArray(v) ? v[0] : v;
        if (typeof currentFieldName === "string") {
            const prefix = subgroupKeyFromFieldName(currentFieldName);
            const matchingLink = links.find((entry) => {
                if (!prefix) {
                    return false;
                }
                return Boolean(entry.link.subgroups[prefix]);
            });
            return matchingLink?.link.name || links[0].link.name;
        }
        return links[0].link.name;
    });

    const handleLinkChange = useCallback((name: string) => {
        setActiveLinkName(name);
    }, []);

    if (links.length === 1) {
        return <LinkToColumnComponentInit {...props} link={links[0]} />;
    }
    return (
        <div>
            <TabHeader
                activeTab={activeLinkName}
                setActiveTab={handleLinkChange}
                tabs={links.map((l) => l.link.name)}
                indicatorTab={committedLinkIndicator}
            />

            {links.map((link) =>
                link.link.name === activeLinkName ? (
                    <LinkToColumnComponentInit
                        key={link.link.name}
                        {...props}
                        link={link}
                    />
                ) : null
            )}
        </div>
    );
});

export default LinkMulti;
