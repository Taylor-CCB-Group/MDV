import { useEffect, useMemo, useState } from "react";
import type DataStore from "../../datastore/DataStore.js";
import TagModel, { type TagColumn } from "../../table/TagModel";
import { BaseDialog } from "../../utilities/Dialog.js";
import { createMdvPortal } from "@/react/react_utils";
import { X } from "lucide-react";
import { observer } from "mobx-react-lite";
import { Autocomplete, TextField } from "@mui/material";


function MTag({tag, tagModel}: {tag: string, tagModel: TagModel}) {
    return (
    <div className="bg-blue-600 p-2 m-1 rounded-lg flex justify-between items-center">
        {tag}
        <button type="button" onClick={() => tagModel.setTag(tag, false)}>
            <X />
        </button>
    </div>
    )
}

function TagInput({tagModel}: {tagModel: TagModel}) {

    return (
        <input className="p-2" placeholder="annotation..." onKeyDown={e => {
            if (e.key === "Enter") {
                if (e.currentTarget.value === '') return;
                tagModel.setTag(e.currentTarget.value, true);
                e.currentTarget.value = '';
            }
        }} />
    )
}


function TagView({dataStore, columnName}: {dataStore: DataStore, columnName: string}) {
    const tagModel = useMemo(() => new TagModel(dataStore, columnName), [dataStore, columnName]);
    const [tagList, setTagList] = useState(tagModel.getTags());
    const [tagsInSelection, setTagsInSelection] = useState(tagModel.getTagsInSelection());
    // we could also distinguish tags that every row in selection has vs only some rows.
    const tagsNotInSelection = [...tagList].filter(tag => !tagsInSelection.has(tag));
    useEffect(() => {
        const listener = tagModel.addListener(() => {
            console.log('tagModel changed');
            setTagList(tagModel.getTags());
            setTagsInSelection(tagModel.getTagsInSelection());
        });
        return () => {
            // probably not relevant now - expect old tagModel to be garbage collected.
            // any remaining references to it are a bug.
            tagModel.removeListener(listener);
        }
    }, [tagModel]);
    
    return (
        <>
        {/* <ChooseColumn column={column} setColumn={setColumn} dataStore={dataStore} /> */}
        <div className="flex flex-col m-4 gap-2">
            <TagInput tagModel={tagModel} />
            {/* <h2>Tags in selection:</h2> */}
            <div className="flex flex-row">
                {[...tagsInSelection].map(tag => <MTag key={tag} tagModel={tagModel} tag={tag}/>)}
                <span className="saturate-0">
                    {[...tagsNotInSelection].map(tag => <MTag key={tag} tagModel={tagModel} tag={tag}/>)}
                </span>
            </div>
            {/* <h2>All tags:</h2>
            <div className="flex flex-row">
                {[...tagList].map(tag => <MTag key={`b${tag}`} tagModel={tagModel} tag={tag}/>)}
            </div> */}
        </div>
        </>
    )
}


const AnnotationDialogComponent = observer(({dataStore}: {dataStore: DataStore}) => {
    const columns = useMemo(() => dataStore.columns.filter(c => c.datatype === 'multitext').map(c => c.name), [dataStore]);
    const [selectedColumn, setSelectedColumn] = useState<string>();
    return (
    <>
        <Autocomplete
        freeSolo
        options={columns}
        onChange={(_, value) => {
            setSelectedColumn(value);
        }}
        renderInput={(params) => (
            <TextField {...params}
            label="Annotation column"
            placeholder="Annotation column name"
            />
        )}
        />
        {selectedColumn && <TagView dataStore={dataStore} columnName={selectedColumn} />}
    </>);
});

class AnnotationDialogReact extends BaseDialog {
    // tagModel: TagModel; //prefer to keep this state in react... but we do need to know what the dataStore is.
    // tagColumn: DataColumn<'text'>;
    // dataModel: DataModel;
    // tagListElement: HTMLDivElement;
    // tagInput: any;
    root: ReturnType<typeof createMdvPortal>;
    constructor(dataStore: DataStore) {
        super({
            title: `Annotate '${dataStore.name}'`,
            width: 400,
            height: 200,
        }, null);
        // this.outer.classList.add('annotationDialog');
        // this.tagModel = new TagModel(dataStore);
        // this.tagModel = tagModel;
        this.root = createMdvPortal(<AnnotationDialogComponent dataStore={dataStore} />, this.dialog, this);
    }
    close(): void {
        super.close();
        this.root.unmount();
    }
}
// https://github.com/Taylor-CCB-Group/MDV/discussions/44
BaseDialog.experiment['AnnotationDialogReact'] = AnnotationDialogReact;
export default 'AnnotationDialogReact loaded';