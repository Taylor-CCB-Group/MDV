import type { ColumnSelectionProps, CTypes } from "@/lib/columnTypeHelpers";
import { observer } from "mobx-react-lite";
import type { DropdownMappedValues, GuiSpec } from "@/charts/charts";
import { useCallback, useEffect, useMemo, useState } from "react";
import { action, makeAutoObservable, observe } from "mobx";
import { DropdownAutocompleteComponent } from "./SettingsDialogComponent";
import { g, isArray } from "@/lib/utils";
import { getFieldName } from "@/links/link_utils";
import { useRowsAsColumnsLinks } from "../chartLinkHooks";

type LocalProps<T extends CTypes, M extends boolean> = ColumnSelectionProps<
    T,
    M
> & {
    link: ReturnType<typeof useRowsAsColumnsLinks>[0];
};


/**
 * A hook for managing the state of manually chosen columns from a rows-as-columns link.
 * 
 * Should encapsulate
 * - making sure that an appropriate type of value is selected
 * - providing a way to select the value that has a simple interface, managing the complexity of incoming generic types
 * 
 * In the end... just bundling everything about state management for this component into a hook, returning a spec.
 */
function useLinkSpec<T extends CTypes, M extends boolean>(props: LocalProps<T, M>) {
    const { linkedDs, link } = props.link;
    const { name_column, name, subgroups } = link;
    const targetColumn = linkedDs.dataStore.columnIndex[name_column];
    if (!targetColumn) {
        throw new Error(`Could not find name_column '${name_column}' in linked dataset`);
    }
    if (linkedDs.size > 2**16) {
        // If we have more than 2^16 rows, I think that means there must be some duplicates...
        // so this is potentially somewhat undefined behaviour.
        console.warn("Linked dataset has more than 2^16 rows, this may lead to unexpected behaviour");
    }
    // we need to associate each value with the row index that references it...
    // as far as the code is concerned, this isn't necessarily a 1:1 mapping - but as my understanding of the logic
    // it should be...
    // We may have
    // - more than one row with the same "gene_id"; I think this leads to undefined behaviour, we should give a warning
    // - "gene_id" that doesn't have any corresponding row... maybe this doesn't matter, but we should filter those values out???
    //   (and probably also log a warning)

    // Loop through the values and create the format in pipes '|', by making use of the index.
    //todo - we need to let the user select subgroup...
    //^^ do we allow a mix-match of links/subgroups? We could... but that doesn't mean we should.
    //(we should only do so if we have a clean way of presenting it to the user)
    const sg = Object.keys(subgroups)[0];
    
    // return values in a format that can be used by the dropdown component
    const values = useMemo(() => {
        // it was possible to get to here with undefined link.valueToRowIndex...
        // should now be avoided by only rendering the component that uses this hook when the link is ready.
        if (!link.valueToRowIndex) {
            // because in certain instances the link isn't ready yet...
            // we can handle that here, but ideally we should do it earlier in the process
            console.warn("Link not ready yet - waiting to get values for UI... should be handled differently.");
            link.initPromise.then(() => {
                console.log(link.valueToRowIndex);
            })
            return [[], "label", "fieldName"] as DropdownMappedValues<"label", "fieldName">;
        }
        // fix issue with safari Map.entries() not being iterable/not having a map method
        const labelValueObjects = [...link.valueToRowIndex.entries()].map(([value, rowIndex]) => ({
            label: value,
            fieldName: getFieldName(sg, value, rowIndex)
        }));
        // tricky type here... but 'satisfies' gives us a bit of confidence that we're getting something that will work
        // ! better than `as` typecasting when possible.
        return [[...labelValueObjects], "label", "fieldName"] satisfies DropdownMappedValues<"label", "fieldName">;
    }, [sg, link.valueToRowIndex, link.initPromise]);


    const { setSelectedColumn } = props;
    const isMultiType = props.multiple;

    // setting internal state will be in the dropdown component...
    // this will be a side-effect to set the state used by the wider application
    const setValue = useCallback((v: string | string[]) => {
        //@ts-ignore could check to the nth degree with props.multiple
        setSelectedColumn(v);
    }, [setSelectedColumn]);
    
    const getSafeInternalValue = useCallback((v: typeof props.current_value) => {
        const defaultVal = props.multiple ? [] : "";
        if (!v) return defaultVal;
        if (props.multiple) {
            // does the current value look like ours? we're not assuming it should:
            // we might not be the active mode...
            if (!isArray(v)) {
                return defaultVal;
            }
            if (v.length === 0) {
                return defaultVal;
            }
            const [firstValue] = v;
            if (typeof firstValue !== 'string') {
                return defaultVal;
            }
            // extracting value
            if (!firstValue.includes("|")) {
                return defaultVal;
            }
            // check that the current value is in the list of possible values
            //! this could be somewhere as a helper function... 
            // given a field name, figure out the row value that would have been used to generate it.
            const val = firstValue.split("|")[1].split(`(${sg})`)[0].trimEnd();
            if (!link.valueToRowIndex?.has(val)) {
                return defaultVal;
            }
            return v;
        } else {
            // if it's an array, we're not the active mode
            if (isArray(v)) {
                return defaultVal;
            }
            if (typeof v !== 'string') {
                return defaultVal;
            }
            if (!v.includes("|")) {
                return defaultVal;
            }
            // extracting value
            const val = v.split("|")[1].split(`(${sg})`)[0].trimEnd();
            // check that the current value is in the list of possible values
            if (!link.valueToRowIndex?.has(val)) {
                return defaultVal;
            }
            return v;
        }
    }, [link.valueToRowIndex, props.multiple, sg]);
    //! there is a bug - when values from a different link are selected, the spec doesn't update
    // we want to be able to set more different combinations of values - which will mean more changes to state management
    // for now, the UI should better reflect what the user is allowed to select.
    const [initialValue] = useState(() => getSafeInternalValue(props.current_value));

    /// we don't want a new spec every time the value changes... we want to give it an observable value
    // (I think)- definition of this spec may want to be in the hook, with an observable value that isn't
    // the real props.current_value???
    const [spec] = useState(() => makeAutoObservable(g<"multidropdown" | "dropdown">({
        type: isMultiType ? 'multidropdown' : 'dropdown',
        label: `specific '${name}'`, //don't really want to have this label here, would take a bigger change to remove
        // no ts error! praise be!
        values,
        current_value: initialValue,
        func: action((v) => {
            // something about having setValue as a dependency makes it end up re-running this...
            // & we really want the spec to be stable.
            //! for some reason if we call this in here we get glitchy behaviour, not really sure
            // setValue(v);
        })
    // })), [values, name, isMultiType, initialValue]);
    })));

    useEffect(() => {
        // do you swim with the mobx stream? or against it?
        // took me a while to realise that this is the direction we should be going in...
        return observe(spec, "current_value", (change) => {
            if (change.type === 'update') {
                setValue(change.newValue);
            }
        });
    }, [spec, setValue]);
    
    return spec;
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
    // is all of the complexity being in the hook really making this component simpler?
    // at least it feels simple if you ignore everything in the hook...
    const spec = useLinkSpec(props);
    return (
        <div className="flex flex-col p-1">
            {/* don't want to have this typecast here, slight nuisance. */}
            <DropdownAutocompleteComponent props={spec as GuiSpec<"dropdown"> | GuiSpec<"multidropdown">} />
        </div>
    )
})

/**
 * Avoid dealing with internal state of the link until it's ready.
 */
const LinkToColumnComponentInit = observer(<T extends CTypes, M extends boolean>(props: LocalProps<T, M>) => {
    const { link } = useRowsAsColumnsLinks()[0];
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

const LinkMulti = observer(
    <T extends CTypes, M extends boolean>(
        props: ColumnSelectionProps<T, M>,
    ) => {
        const links = useRowsAsColumnsLinks();
        const [activeLinkIndex, setActiveLinkIndex] = useState(() => {
            // Determine the initial active link based on `current_state` and `FieldName`
            // this logic might be adjusted if we allow combined selections with multiple links
            // also we might want to persist as user changes tab outside of this component
            const v = props.current_value;
            const currentFieldName = isArray(v) ? v[0] : v;
            if (typeof currentFieldName === "string") {
                const matchingIndex = links.findIndex((link) => {
                    const sgName = Object.keys(link.link.subgroups)[0];
                    return currentFieldName.startsWith(`${sgName}|`);
                });
                return matchingIndex !== -1 ? matchingIndex : 0;
            }
            return 0;
        });

        const handleLinkChange = (index: number) => {
            setActiveLinkIndex(index);
        };
        // don't show tab ui if we only have one link
        if (links.length === 1) {
            return <LinkToColumnComponentInit {...props} link={links[0]} />;
        }
        return (
            <div>
                {/* Render tabs or radio buttons for selecting the active link */}
                <div className="flex space-x-2 mb-4">
                    {links.map((link, index) => (
                        <button
                            key={link.link.name}
                            className={`px-4 py-1 border-b-2 ${
                                index === activeLinkIndex
                                    ? "font-bold"
                                    : "font-light"
                            }`}
                            onClick={() => handleLinkChange(index)}
                        >
                            {link.link.name}
                        </button>
                    ))}
                </div>

                {/* Render the active link's UI */}
                {links.map((link, index) =>
                    index === activeLinkIndex ? (
                        <LinkToColumnComponentInit
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


export default LinkMulti;