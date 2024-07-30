import { useEffect, useState } from "react";
import type DataStore from "../../datastore/DataStore.js";
import type TagModel from "../../table/TagModel";
import { BaseDialog } from "../../utilities/Dialog.js";
import { createMdvPortal } from "@/react/react_utils";
import { X } from "lucide-react";


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


function AnnotationDialogComponent({tagModel}: {tagModel: TagModel}) {
    const [tagList, setTagList] = useState(tagModel.getTags());
    const [tagsInSelection, setTagsInSelection] = useState(tagModel.getTagsInSelection());
    useEffect(() => {
        const listener = tagModel.addListener(() => {
            console.log('tagModel changed');
            setTagList(tagModel.getTags());
            setTagsInSelection(tagModel.getTagsInSelection());
        });
        return () => {
            tagModel.removeListener(listener);
        }
    })
    
    return (
        <>
        <div className="flex flex-col m-4 gap-2">
            <TagInput tagModel={tagModel} />
            <h2>Tags in selection:</h2>
            <div className="flex flex-row">
                {[...tagsInSelection].map(tag => <MTag key={`a${tag}`} tagModel={tagModel} tag={tag}/>)}
            </div>
            {/* <h2>All tags:</h2>
            <div className="flex flex-row">
                {[...tagList].map(tag => <MTag key={'b'+tag} tagModel={tagModel} tag={tag}/>)}
            </div> */}
        </div>
        </>
    )
}


class AnnotationDialogReact extends BaseDialog {
    tagModel: TagModel;
    // tagColumn: DataColumn<'text'>;
    // dataModel: DataModel;
    // tagListElement: HTMLDivElement;
    // tagInput: any;
    root: ReturnType<typeof createMdvPortal>;
    constructor(dataStore: DataStore, tagModel: TagModel) {
        super({
            title: "Annotate selection",
            width: 400,
            height: 200,
        }, null);
        this.outer.classList.add('annotationDialog');
        // this.tagModel = new TagModel(dataStore);
        this.tagModel = tagModel;
        this.root = createMdvPortal(<AnnotationDialogComponent tagModel={this.tagModel} />, this.dialog);
    }
    close(): void {
        super.close();
        this.root.unmount();
    }
}
// https://github.com/Taylor-CCB-Group/MDV/discussions/44
BaseDialog.experiment['AnnotationDialogReact'] = AnnotationDialogReact;
export default 'AnnotationDialogReact loaded';