import { useEffect, useState } from "react";
import DataStore from "../../datastore/DataStore.js";
import TagModel from "../../table/TagModel";
import { BaseDialog } from "../../utilities/Dialog.js";
import { createMdvPortal } from "@/react/react_utils";

function AnnotationDialogComponent({tagModel}: {tagModel: TagModel}) {
    const [tagList, setTagList] = useState(tagModel.getTags());
    const [tagsInSelection, setTagsInSelection] = useState(tagModel.getTagsInSelection());
    useEffect(() => {
        const listener = tagModel.addListener(() => {
            setTagList(tagModel.getTags());
            setTagsInSelection(tagModel.getTagsInSelection());
        });
        return () => {
            tagModel.removeListener(listener);
        }
    })
    
    return (
        <div>
            <input></input>
            <h2>Tags in selection:</h2>
            <ul>
                {[...tagsInSelection].map(tag => <li key={'a'+tag}>{tag}</li>)}
            </ul>
            <h2>All tags:</h2>
            <ul>
                {[...tagList].map(tag => <li key={'b'+tag}>{tag}</li>)}
            </ul>
        </div>
    )
}


class AnnotationDialogReact extends BaseDialog {
    tagModel: TagModel;
    // tagColumn: DataColumn<'text'>;
    // dataModel: DataModel;
    tagListElement: HTMLDivElement;
    tagInput: any;
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