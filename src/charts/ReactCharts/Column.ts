import { types,flow } from "mobx-state-tree";
function getColumn(datastore: string, field: string) {
    return new Promise((resolve, reject) => {
        const cm = window.mdv.chartManager;
        cm._getColumnsThen(datastore, [field],() => {
            resolve(cm.getDataSource(datastore).columnIndex[field]);
        })
                        
    });
}


const LinkModel = types.model({
    datastore: types.string,
    subgroup: types.string,
});

const Column = types
    .model({
        datastore: types.string,
        field: types.string,
        link: types.maybe(LinkModel)
    })
    .actions((self) => ({
        setDataStore(datastore: string) {
            //checks here
            self.datastore = datastore;
        },
        setField(field: string) {
            //checks here
            self.field = field;
        },
        setLink(datastore: string, subgroup: string) {
            //checks here
            self.link = LinkModel.create({ datastore, subgroup });
        },
        removeLink() {
            self.link = undefined;
        }
    }));

export default Column


