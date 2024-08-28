import { useEffect, useMemo, useState } from "react";
import type DataStore from "../../datastore/DataStore.js";
import TagModel, { type TagColumn } from "../../table/TagModel";
import { BaseDialog } from "../../utilities/Dialog.js";
import { createMdvPortal } from "@/react/react_utils";
import { X } from "lucide-react";
import { observer } from "mobx-react-lite";
import TextField from "@mui/material/TextField";
import Autocomplete from "@mui/material/Autocomplete";
import Checkbox from "@mui/material/Checkbox";
import Chip from "@mui/material/Chip";
import CheckBoxOutlineBlankIcon from "@mui/icons-material/CheckBoxOutlineBlank";
import CheckBoxIcon from "@mui/icons-material/CheckBox";

function MTag({ tag, tagModel }: { tag: string; tagModel: TagModel }) {
    return (
        <div className="bg-blue-600 p-2 m-1 rounded-lg flex justify-between items-center">
            {tag}
            <button type="button" onClick={() => tagModel.setTag(tag, false)}>
                <X />
            </button>
        </div>
    );
}

function TagInput({ tagModel }: { tagModel: TagModel }) {
    return (
        <input
            className="p-2"
            placeholder="annotation..."
            onKeyDown={(e) => {
                if (e.key === "Enter") {
                    if (e.currentTarget.value === "") return;
                    tagModel.setTag(e.currentTarget.value, true);
                    e.currentTarget.value = "";
                }
            }}
        />
    );
}

function TagViewX({
    dataStore,
    columnName,
}: { dataStore: DataStore; columnName: string }) {
    const tagModel = useMemo(
        () => new TagModel(dataStore, columnName),
        [dataStore, columnName],
    );
    const [tagList, setTagList] = useState(tagModel.getTags());
    const [tagsInSelection, setTagsInSelection] = useState(
        tagModel.getTagsInSelection(),
    );
    // we could also distinguish tags that every row in selection has vs only some rows.
    const tagsNotInSelection = [...tagList].filter(
        (tag) => !tagsInSelection.has(tag),
    );
    useEffect(() => {
        const listener = tagModel.addListener(() => {
            console.log("tagModel changed");
            setTagList(tagModel.getTags());
            setTagsInSelection(tagModel.getTagsInSelection());
        });
        return () => {
            // probably not relevant now - expect old tagModel to be garbage collected.
            // any remaining references to it are a bug.
            tagModel.removeListener(listener);
        };
    }, [tagModel]);

    return (
        <>
            {/* <ChooseColumn column={column} setColumn={setColumn} dataStore={dataStore} /> */}
            <div className="flex flex-col m-4 gap-2">
                <TagInput tagModel={tagModel} />
                {/* <h2>Tags in selection:</h2> */}
                <div className="flex flex-row">
                    {[...tagsInSelection].map((tag) => (
                        <MTag key={tag} tagModel={tagModel} tag={tag} />
                    ))}
                    <span className="saturate-0">
                        {[...tagsNotInSelection].map((tag) => (
                            <MTag key={tag} tagModel={tagModel} tag={tag} />
                        ))}
                    </span>
                </div>
                {/* <h2>All tags:</h2>
            <div className="flex flex-row">
                {[...tagList].map(tag => <MTag key={`b${tag}`} tagModel={tagModel} tag={tag}/>)}
            </div> */}
            </div>
        </>
    );
}

const icon = <CheckBoxOutlineBlankIcon fontSize="small" />;
const checkedIcon = <CheckBoxIcon fontSize="small" />;

function TagView({
    dataStore,
    columnName,
}: { dataStore: DataStore; columnName: string }) {
    const tagModel = useMemo(
        () => new TagModel(dataStore, columnName),
        [dataStore, columnName],
    );
    useEffect(() => {
        const listener = tagModel.addListener(() => {
            console.log("tagModel changed");
            setTagList(tagModel.getTags());
            setTagsInSelection(tagModel.getTagsInSelection());
            console.log("tagList", tagModel.getTags());
            console.log("tagsInSelection", tagModel.getTagsInSelection());
        });
        return () => {
            // probably not relevant now - expect old tagModel to be garbage collected.
            // any remaining references to it are a bug.
            // ^^^ also the design is likely to change and there probably are bugs ATM ^^^
            tagModel.removeListener(listener);
        };
    }, [tagModel]);

    const [tagList, setTagList] = useState(
        tagModel.isReady ? tagModel.getTags() : new Set<string>(),
    );
    const [tagsInSelection, setTagsInSelection] = useState(
        tagModel.isReady ? tagModel.getTagsInSelection() : new Set<string>(),
    );

    return (
        <>
            <h3>
                Annotations (on {tagModel.dataModel.getLength()} selected rows):
            </h3>
            <Autocomplete
                freeSolo
                multiple
                options={[...tagList]} //strings vs objects? strings work more easily with freeSolo
                value={[...tagsInSelection]}
                onChange={(_, value) => {
                    // we only handle items that have been either removed or added here
                    // - items that were previously only present in part of the selection
                    // but have been toggled to 'on' will be handled elsewhere (onClick in renderOption components)
                    const addedTags = value.filter(
                        (v) => !tagsInSelection.has(v),
                    );
                    for (const t of addedTags) {
                        console.log("main Autocomplete onChange adding", t);
                        tagModel.setTag(t, true);
                    }
                    const removedTags = [...tagsInSelection].filter(
                        (v) => !value.includes(v),
                    );
                    for (const t of removedTags) {
                        console.log("main Autocomplete onChange removing", t);
                        tagModel.setTag(t, false);
                    }
                }}
                renderInput={(props) => {
                    const { key, ...p } = props as typeof props & {
                        key: string;
                    }; //questionable mui types?
                    return (
                        <TextField
                            key={key}
                            {...p}
                            label={`add tags to '${columnName}'...`}
                            onKeyDown={(e) => {
                                if (e.key === "Enter") {
                                    const { target } = e;
                                    if (target instanceof HTMLInputElement) {
                                        tagModel.setTag(target.value, true);
                                    }
                                }
                            }}
                        />
                    );
                }}
                // getOptionLabel={o => typeof o === "string" ? o : o.tag}
                // getOptionLabel={o => o}
                renderOption={(props, tag, { selected }) => {
                    const { key, ...optionProps } = props as typeof props & {
                        key: string;
                    }; //questionable mui types?
                    const indeterminate =
                        !tagModel.entireSelectionHasTag(tag) &&
                        tagsInSelection.has(tag);
                    return (
                        <li
                            key={key}
                            {...optionProps}
                            onClick={() => {
                                if (indeterminate) tagModel.setTag(tag, true);
                                else tagModel.setTag(tag, !selected);
                            }}
                        >
                            <Checkbox
                                icon={icon}
                                checkedIcon={checkedIcon}
                                style={{ marginRight: 8 }}
                                // checked={option.inSelection === 'entire'}
                                // indeterminate={option.inSelection === 'partial'}
                                checked={tagModel.entireSelectionHasTag(tag)}
                                indeterminate={indeterminate}
                            />
                            {tag}
                        </li>
                    );
                }}
                renderTags={(tagValue, getTagValues) => {
                    //seems to be a material-ui bug with not properly handling key / props...
                    //https://stackoverflow.com/questions/75818761/material-ui-autocomplete-warning-a-props-object-containing-a-key-prop-is-be
                    return tagValue.map((option, index) => {
                        const { key, ...tagValues } = getTagValues({ index });
                        return <Chip key={key} {...tagValues} label={option} />;
                    });
                }}
            />
        </>
    );
}

const AnnotationDialogComponent = observer(
    ({ dataStore }: { dataStore: DataStore }) => {
        const columns = useMemo(
            () =>
                dataStore.columns
                    .filter((c) => c.datatype === "multitext")
                    .map((c) => c.name),
            [dataStore],
        );
        const [selectedColumn, setSelectedColumn] = useState<string>();
        type NameValid = "ok" | "empty" | "clash";
        const [isValidName, setIsValidName] = useState<NameValid>("empty");
        const validateInput = (value: string) => {
            if (value.trim() === "") {
                setIsValidName("empty");
                return false;
            }
            const existingCol = dataStore.columnIndex[value];
            if (existingCol && existingCol.datatype !== "multitext") {
                setIsValidName("clash");
                return false;
            }
            setIsValidName("ok");
            return true;
        };
        return (
            <>
                <Autocomplete
                    freeSolo
                    options={columns}
                    onChange={(_, value) => {
                        if (validateInput(value)) setSelectedColumn(value);
                    }}
                    onInputChange={({ currentTarget }) => {
                        if (currentTarget instanceof HTMLInputElement) {
                            const { value } = currentTarget;
                            validateInput(value);
                        }
                    }}
                    renderInput={(params) => {
                        const { key, ...p } = params as typeof params & {
                            key: string;
                        };
                        return (
                            <TextField
                                key={key}
                                {...p}
                                color={isValidName ? undefined : "warning"}
                                label="Annotation column"
                                placeholder="Annotation column name"
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
                {isValidName === "clash" && (
                    <p>Incompatible column with that name already exists...</p>
                )}
                {selectedColumn && (
                    <TagView
                        dataStore={dataStore}
                        columnName={selectedColumn}
                    />
                )}
                {!selectedColumn && (
                    <h2>
                        Select or add a column above to use for annotations.
                    </h2>
                )}
            </>
        );
    },
);

class AnnotationDialogReact extends BaseDialog {
    // tagModel: TagModel; //prefer to keep this state in react... but we do need to know what the dataStore is.
    // tagColumn: DataColumn<'text'>;
    // dataModel: DataModel;
    // tagListElement: HTMLDivElement;
    // tagInput: any;
    root: ReturnType<typeof createMdvPortal>;
    constructor(dataStore: DataStore) {
        super(
            {
                title: `Annotate '${dataStore.name}'`,
                width: 400,
                height: 300,
            },
            null,
        );
        // this.outer.classList.add('annotationDialog');
        // this.tagModel = new TagModel(dataStore);
        // this.tagModel = tagModel;
        this.root = createMdvPortal(
            <AnnotationDialogComponent dataStore={dataStore} />,
            this.dialog,
            this,
        );
    }
    close(): void {
        super.close();
        this.root.unmount();
    }
}
// https://github.com/Taylor-CCB-Group/MDV/discussions/44
BaseDialog.experiment["AnnotationDialogReact"] = AnnotationDialogReact;
export default "AnnotationDialogReact loaded";
