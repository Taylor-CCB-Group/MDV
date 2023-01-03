import DataStore from "../../datastore/DataStore.js";
import TagModel from "../../table/TagModel";
import {BaseDialog} from "../../utilities/Dialog.js";
import {createEl} from "../../utilities/Elements.js";


export default class AnnotationDialog extends BaseDialog {
    tagModel: TagModel;
    // tagColumn: DataColumn<'text'>;
    // dataModel: DataModel;
    tagListElement: HTMLDivElement;
    tagInput: any;
    constructor(dataStore: DataStore) {
        super({
            title: "Annotate selection",
            width: 400,
            height: 200,
            columns: 1,
            buttons: [
                {
                    text: "Add tag to selection",
                    method: () => {
                        const tag = this.tagInput.value;
                        this.setTag(tag, true);
                        this.updateTagList();
                    }
                },
                {
                    text: "Remove tag from selection",
                    method: () => {
                        const tag = this.tagInput.value;
                        this.setTag(tag, false);
                        this.updateTagList();
                    }
                },
            ]
        }, null);
        this.outer.classList.add('annotationDialog');
        this.tagModel = new TagModel(dataStore);
        this.tagModel.addListener(() => this.updateTagList());
        

        const tagInput = createEl("input", {}, this.columns[0]); //chrome treats as password field (but without hidden value)?
        tagInput.addEventListener('keydown', e => {
            if (e.key == "Enter") this.setTag(tagInput.value, true);
        });
        tagInput.focus();
        this.tagInput = tagInput;
        // const addTagButton = createEl("button", {text: "Add tag to selection", classes:["ciview-button"]}, this.columns[1]);
        // addTagButton.addEventListener("click", () => {
        //     setTagOnAllSelectedValues(tagInput.value, this.tagColumn, dataModel);
        //     this.updateTagList();
        // });
        this.tagListElement = createEl("div", {}, this.columns[0]);
        this.updateTagList();
    }
    setTag(tag: string, tagValue: boolean) {
        this.tagModel.setTag(tag, tagValue);
        this.updateTagList();
    }
    updateTagList() {
        const {tagModel} = this;
        const parent = this.tagListElement;
        parent.innerHTML = '';
        const allTagsEl = createEl('ul', {}, parent);
        for (const tag of tagModel.getTags()) {
            const el = createEl('li', {}, allTagsEl);
            const button = createEl('button', {text: tag}, el); //more accessible... why even have the ul/li?
            button.addEventListener('click', () => this.setTag(tag, true));
        }
        const selTagsEl = createEl('ul', {}, parent);
        for (const tag of tagModel.getTagsInSelection()) {
            const all = this.tagModel.entireSelectionHasTag(tag);
            const li = createEl('li', {classes: [all ? "allTagged" : "someTagged"]}, selTagsEl);
            const button = createEl('button', {text: tag}, li);
            //bug: this is removing other tags, not just this one
            button.addEventListener('click', () => this.setTag(tag, !all));
        }
    }
}
