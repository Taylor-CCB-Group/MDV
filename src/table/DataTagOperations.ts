import { DataColumn } from "../charts/charts";

/** Treating `col.values` as strings containing space-separated 'tags',
 * find all the indices that include 'tag' as one of those tags.
 */
export function getTagValueIndices(tag: string, col: DataColumn<'text'>) {
    // if we're going to treat values this way, we might want to consider escaping spaces,
    // or treating the value passed to 'tag' as potentially containing multiple tags?
    // we could use some horrible special character as tag delimiter
    if (tag.split(' ').length > 1) console.error('todo: process space in tag input');
    return col.values.map((v, i) => v.split(' ').includes(tag) ? i : -1).filter(i => i != -1);
}

type DataModel = {data: Int32Array, _getValueIndex(value: string, col: DataColumn<'text'>)}

export function setTagOnAllValues(tag: string, col: DataColumn<'text'>, dataModel: DataModel) {
    const indicesWithTag = getTagValueIndices(tag, col); //refers to values that already contain 'tag'

    //map from index of value without tag, to index of that value with the tag added, if it already exists...
    //or -1 if it doesn't (yet) exist.
    //As we process the list, every time we have a new combination of tags (a new value), 
    //the mapping will be added to
    const mapToAppendedTag = new Map<number, number>();
    col.values.filter((_, i) => !indicesWithTag.includes(i)).forEach((v, i) => {
        const tags = [tag, ...v.split(' ')].sort(); //sort *after* adding tag
        //we could assert !tags.includes(tag) before it was added
        const valueWithTag = tags.join(' ');
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
            const taggedValue = [tag, ...untaggedValue.split(' ')].sort().join(' ');
            taggedIndex = dataModel._getValueIndex(taggedValue, col);
            mapToAppendedTag.set(untaggedIndex, taggedIndex);
        }
        col.data[data[i]] = taggedIndex;
    }
    console.log(`updated ${count}/${data.length} selected rows`);
}
