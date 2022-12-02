import { DataColumn } from "../charts/charts";

type TagColumn = DataColumn<'text'>;

const SEP = /\W*\,\W*/; //separate by comma with whitespace trimmed
const JOIN = ', '; //join with comma and space.

const splitTags = (value: string) => value.split(SEP).filter(v=>v);


/** Treating `col.values` as strings containing comma-separated 'tags',
 * find all the indices that include 'tag' as one of those tags.
 */
export function getTagValueIndices(tag: string, col: TagColumn) {
    if (tag.split(SEP).length > 1) console.error('todo: process multiple tags in tag input');
    return col.values.map((v, i) => v.split(SEP).includes(tag) ? i : -1).filter(i => i != -1);
}

type DataModel = {data: Int32Array, _getValueIndex(value: string, col: TagColumn)}

function getValueIndex(value: string, col: TagColumn) {
    let i = col.values.indexOf(value);
    if (i === -1) {
        col.values.push(value);
        i = col.values.length - 1;
    }
    return i;
}

export function setTagOnAllValues(tag: string, col: TagColumn, dataModel: DataModel) {
    sanitizeTags(col);
    const indicesWithTag = getTagValueIndices(tag, col); //refers to values that already contain 'tag'

    //map from index of value without tag, to index of that value with the tag added, if it already exists...
    //or -1 if it doesn't (yet) exist.
    //As we process the list, every time we have a new combination of tags (a new value), 
    //the mapping will be added to
    const mapToAppendedTag = new Map<number, number>();
    col.values.filter((_, i) => !indicesWithTag.includes(i)).forEach((v, i) => {
        const tags = [tag, ...splitTags(v)].sort(); //sort *after* adding tag
        //we could assert !tags.includes(tag) before it was added
        const valueWithTag = tags.join(JOIN);
        mapToAppendedTag.set(i, col.values.indexOf(valueWithTag));
    });

    const {data} = dataModel;
    let count = 0;
    for (let i=0; i<data.length; i++) {
        // if data was a bitmap, we'd just "or" our value
        // so things like searching, adding tag would be very quick, but we might use a bit more space
        // maybe we'd probably need more changes... but perhaps it would be clearer separation
        // to introduce a new datatype with explicit semantics than hacking on top of 'text'?
        const currentVal = col.data[data[i]];
        
        if (indicesWithTag.includes(currentVal)) continue;
        count++;
        // we've found a row that doesn't have the tag... if we were to add the tag, 
        // would that be a new value (ie, nothing else would have that combination of tags)?
        const untaggedIndex = currentVal;
        let taggedIndex = mapToAppendedTag.get(untaggedIndex);
        if (taggedIndex == -1) {
            const untaggedValue = col.values[untaggedIndex];
            const taggedValue = [tag, ...splitTags(untaggedValue)].sort().join(JOIN);
            taggedIndex = dataModel._getValueIndex(taggedValue, col) as number;
            mapToAppendedTag.set(untaggedIndex, taggedIndex);
        }
        col.data[data[i]] = taggedIndex!;
    }
    console.log(`updated ${count}/${data.length} selected rows`);
}

export function removeTagFromAllSelectedValues(tag: string, col: TagColumn, dataModel: DataModel) {
    // it may be possible that certain tags previously only appeared in combination with others,
    // so we should potentially cleanup... sanitizeTags doesn't currently remove unreferenced values.
}

/** processes a given column so that tags appear in sorted order and without repitition.
 * Note that this may potentially alter data in rows that where the corresponding value already
 * satisfied these conditions.
 */
export function sanitizeTags(col: TagColumn) {
    const vals = col.values.map((unsorted, i) => {
        const sorted = [...new Set(splitTags(unsorted))].sort().join(JOIN);
        console.log({unsorted, sorted});
        return {i, unsorted, sorted}
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
}

export function getTags(col: TagColumn) {
    return new Set(col.values.map(splitTags).flat());
}

export function getTagsInSelection(col: TagColumn, dataModel: DataModel) {
    const usedValuesByIndex = new Set<number>();
    dataModel.data.forEach(i => usedValuesByIndex.add(col.data[i]));
    const values = [...usedValuesByIndex].map(i => col.values[i]);
    return new Set(values.map(splitTags).flat());
}