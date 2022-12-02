import DataStore from "../../datastore/DataStore.js";
import TagModel from "../../table/TagModel";
import {BaseDialog} from "../../utilities/Dialog.js";
import {createEl} from "../../utilities/Elements.js";


export default class AnnotationDialog extends BaseDialog {
    tagModel: TagModel;
    // tagColumn: DataColumn<'text'>;
    // dataModel: DataModel;
    tagList: HTMLDivElement;
    tagInput: any;
    constructor(dataStore: DataStore) {
        super({
            title: "Annotate selection",
            width: 400,
            height: 400,
            columns: 2,
            buttons: [
                {
                    text: "Add tag to selection",
                    method: () => {
                        const tag = this.tagInput.value;
                        this.tagModel.setTag(tag, true);
                        this.updateTagList();
                    }
                },
                {
                    text: "Remove tag from selection",
                    method: () => {
                        const tag = this.tagInput.value;
                        this.tagModel.setTag(tag, false);
                        this.updateTagList();
                    }
                },
            ]
        }, null);
        this.outer.classList.add('annotationDialog');
        this.tagModel = new TagModel(dataStore);
        this.tagModel.addListener(() => this.updateTagList());
        

        const tagInput = createEl("input", {}, this.columns[0]); //chrome treats as password field (but without hidden value)?
        this.tagInput = tagInput;
        // const addTagButton = createEl("button", {text: "Add tag to selection", classes:["ciview-button"]}, this.columns[1]);
        // addTagButton.addEventListener("click", () => {
        //     setTagOnAllSelectedValues(tagInput.value, this.tagColumn, dataModel);
        //     this.updateTagList();
        // });
        this.tagList = createEl("div", {}, this.columns[0]);
        this.updateTagList();
    }
    updateTagList() {
        const allTags = [...this.tagModel.getTags()];
        const els = allTags.map(tag => `<li>${tag}</li>`);
        const renderTag = (tag: string) => {
            const all = this.tagModel.entireSelectionHasTag(tag); //needs fixing.
            const style = all ? `style="border 1px solid white"` : ""
            return `<li ${style}>${tag}</li>`;
        }
        const tm = this.tagModel;
        this.tagList.innerHTML = `<h2>all tags:</h2>
            <ul>${els.join('')}</ul>
            <h2>tags used in selection:</h2>
            <ul>
            ${[...tm.getTagsInSelection()].map(renderTag).join('')}
            </ul>
        `;
    }
}
