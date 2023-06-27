import { DataColumn } from "../charts/charts";
import DataStore from "../datastore/DataStore";
import { DataModel } from "./DataModel";

type TagColumn = DataColumn<'multitext'>; //multitext...?

const SEP = /\W*\;\W*/; //separate by semi-colon with whitespace trimmed
const JOIN = '; '; //join with semi-colon and space.

const splitTags = (value: string) => value.split(SEP).filter(v=>v);


/** Treating `col.values` as strings containing semi-colon-separated 'tags',
 * find all the indices that include 'tag' as one of those tags.
 */
function getTagValueIndices(tag: string, col: TagColumn) {
    if (tag.split(SEP).length > 1) console.error('todo: process multiple tags in tag input');
    return col.values.map((v, i) => v.split(SEP).includes(tag) ? i : -1).filter(i => i != -1);
}

export default class TagModel {
    readonly tagColumn: TagColumn;
    readonly dataModel: DataModel;
    listeners: (()=>void)[] = [];
    constructor(dataStore: DataStore) {
        this.dataModel = new DataModel(dataStore, {autoupdate: true});
        this.dataModel.updateModel();
        //nb, this will replace any other "tag" listener on model (but the model is local to this object)
        this.dataModel.addListener("tag", () => {
            this.listeners.map(callback => callback());
        })
        if (!dataStore.columnIndex['__tags']) this.dataModel.createColumn('__tags', null);
        this.tagColumn = dataStore.columnIndex['__tags'];
    }
    addListener(callback: ()=>void) {
        this.listeners.push(callback);
        return callback;
    }
    removeListener(callback: ()=>void) {
        const i = this.listeners.indexOf(callback);
        if (i != -1) this.listeners.splice(i, 1);
    }
    setTag(tag: string, tagValue = true) {
        this.dataModel.updateModel();
        setTag({tag, ...this}, tagValue);
    }
    getTags() {
        return getTags(this.tagColumn);
    }
    entireSelectionHasTag(tag: string) {
        const col = this.tagColumn;
        const tagIndices = getTagValueIndices(tag, col);
        return !this.dataModel.data.some(v => !tagIndices.includes(col.data[v]));
    }
    getTagsInSelection() {
        const {dataModel, tagColumn: col} = this;
        const usedValuesByIndex = new Set<number>();
        dataModel.data.forEach(i => usedValuesByIndex.add(col.data[i]));
        const values = [...usedValuesByIndex].map(i => col.values[i]);
        return new Set(values.map(splitTags).flat());
    }
}

function getValueIndex(value: string, col: TagColumn) {
    let i = col.values.indexOf(value);
    if (i === -1) {
        col.values.push(value);
        if (col.values.length > 256) throw new Error(`text column '${col.name}' exceeded 256 values when adding '${value}'`);
        i = col.values.length - 1;
    }
    return i;
}
function setTag(obj: {tag: string, tagColumn: TagColumn, dataModel: DataModel}, tagValue = true) {
    const {tag, tagColumn, dataModel} = obj;
    setTagOnAllSelectedValues(tag, tagColumn, dataModel, true, tagValue);
}
function setTagOnAllSelectedValues(tagToChange: string, col: TagColumn, dataModel: DataModel, notify = true, tagValue = true) {
    // if (dataModel.data.length == 0) {
    //     setTagOnAllValues(tag, col, dataModel, notify); //don't want to maintain permutations like this.
    //     return;
    // }
    sanitizeTags(col);
    const indicesWithTag = getTagValueIndices(tagToChange, col); //refers to values that already contain 'tag'

    //If tagValue is true, map from the index of the value without the tag, to the index of the value with the tag.
    //If false, the opposite.
    //-1 if a value with corresponding set of tags doesn't (yet) exist.
    //As we process the list, every time we need to use a new combination of tags (a new value), those -1 values will be set,
    //but we don't add new values until we know we need them.
    const mapToAlteredTag = new Map<number, number>();
    if (tagValue) {
        col.values.filter((_, i) => !indicesWithTag.includes(i)).forEach((v, i) => {
            const tags = [tagToChange, ...splitTags(v)].sort(); //sort *after* adding tag
            //we could assert !tags.includes(tag) before it was added
            const valueWithTag = tags.join(JOIN);
            mapToAlteredTag.set(i, col.values.indexOf(valueWithTag));
        });
    } else {
        col.values.filter((_, i) => indicesWithTag.includes(i)).forEach((v, i) => {
            const tags = splitTags(v).filter(t => t !== tagToChange).sort();
            const valueWithoutTag = tags.join(JOIN);
            mapToAlteredTag.set(i, col.values.indexOf(valueWithoutTag));
        });
    }

    const {data} = dataModel;
    let count = 0;
    for (let i=0; i<data.length; i++) {
        // if data was a bitmap, we'd just do binary operations on our value
        // so things like searching, adding tag would be very quick, but we might use a bit more space
        // maybe we'd probably need more changes... but perhaps it would be clearer separation
        // to introduce a new datatype with explicit semantics than hacking on top of 'text'?
        // this method is a bit unintuitive to write, and could be a source of bugs.
        const currentVal = col.data[data[i]];
        
        if (tagValue == indicesWithTag.includes(currentVal)) continue;
        count++;
        // we've found a row that doesn't have the tag... if we were to add the tag, 
        // would that be a new value (ie, nothing else would have that combination of tags)?
        const untaggedIndex = currentVal;
        let taggedIndex = mapToAlteredTag.get(untaggedIndex);
        // undefined should be never, but it looks like that logic isn't right when removing tags
        if (taggedIndex == undefined || taggedIndex == -1) {
            const untaggedValue = col.values[untaggedIndex];
            const tags = splitTags(untaggedValue);
            const alteredTags = tagValue ? [tagToChange, ...tags] : tags.filter(t => t !== tagToChange);
            const taggedValue = alteredTags.sort().join(JOIN);
            taggedIndex = getValueIndex(taggedValue, col);
            mapToAlteredTag.set(untaggedIndex, taggedIndex);
        }
        col.data[data[i]] = taggedIndex!;
    }
    // console.log(`updated ${count}/${data.length} selected rows`);
    if (notify && count) dataModel.dataStore.dataChanged([col.name]);
}

/** processes a given column so that tags appear in sorted order and without repitition.
 * Note that this may potentially alter data in rows that where the corresponding value already
 * satisfied these conditions.
 */
function sanitizeTags(col: TagColumn, notify = false) {
    const vals = col.values.map((unsorted, i) => {
        const sorted = [...new Set(splitTags(unsorted))].sort().join(JOIN);
        return {i, unsorted, sorted};
    });
    const alreadySorted = vals.filter(v => v.sorted == v.unsorted);
    if (alreadySorted.length === vals.length) return;
    console.log('sanitizing tags...'); //not expected to happen unless other non-tag operations are performed on col?
    const values = alreadySorted.map(v => v.sorted);
    col.values = values;
    
    const mapToSorted = new Map<number, number>();
    // we go over all vals, as values that are to be kept (alreadySorted) may still have their indices changed.
    for (const v of vals) {
        const {i, sorted} = v;
        mapToSorted.set(i, getValueIndex(sorted, col));
    }
    const usedValuesByIndex = new Set<number>();
    for (let i=0; i<col.data.length; i++) {
        const j = col.data[i] = mapToSorted.get(col.data[i])!;
        usedValuesByIndex.add(j);
    }
    const nUnused = vals.length - usedValuesByIndex.size;
    if (nUnused) {
        // ugh, we may want to do another pass... make sure to actually test when I implement this.
        console.warn(`sanitizeTags left ${nUnused} unused values...`);
    }
    if (notify) console.warn('sanitizeTags notify ignored');
}

function getTags(col: TagColumn) {
    return new Set(col.values.map(splitTags).flat());
}

function getTagsInSelection(col: TagColumn, dataModel: DataModel) {
    const usedValuesByIndex = new Set<number>();
    dataModel.data.forEach(i => usedValuesByIndex.add(col.data[i]));
    const values = [...usedValuesByIndex].map(i => col.values[i]);
    return new Set(values.map(splitTags).flat());
}