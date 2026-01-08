import BaseChart from "./BaseChart";
import { createEl } from "../utilities/Elements.js";
import CustomDialog from "./dialogs/CustomDialog.js";
import { useDataSources } from "@/react/hooks";
import {loadColumn} from "@/dataloaders/DataLoaderUtil"


class GenomeBrowser extends BaseChart {
    constructor(dataStore, div, config) {
        super(dataStore, div, config);
        const c = this.config;
        c.view_margins = c.view_margins || { type: "percentage", value: 20 };
        this.contentDiv.style.width = "100%";
        this.contentDiv.style.height = "100%";
        c.type = "genome_browser";
        const add_ruler = !c.tracks.find((x) => x["format"] === "ruler");

        import("../browser/panel.js?v0.3").then(({ default: MLVPanel }) => {
            this.browser = new MLVPanel(
                this.contentDiv,
                {
                    allow_user_interactions: true,
                    show_scale: true,
                    ruler_track: add_ruler,
                },
                c.tracks,
            );

            this.bamscatrack = this.browser.getTrack("_atac_bam_track");
            this.baseTrack = this.browser.getTrack("_base_track");
            this.baseTrack.dataStore = this.dataStore;

            this.addMenuIcon(
                "fa fa-eye-slash",
                config.legend || "Show/Hide Tracks",
                {
                    func: (e) => this.showHideDialog(),
                },
            );
            this.locText = createEl(
                "input",
                { styles: { fontSize: "11px" } },
                this.menuSpace,
            );
            this.locText.addEventListener("keypress", (e) => {
                if (e.key === "Enter") {
                    const loc = this._calculatePosition(this.locText.value);
                    this.browser.update(loc.chr, loc.start, loc.end);
                }
            });
            this.browser.addListener(this.config.id, (type, data, e) =>
                this.onBrowserAction(type, data, e),
            );
            if (c.feature_label) {
                this.setLabelFunction(c.feature_label);
            }
            if (c.color_by) {
                this.colorByColumn(c.color_by, false);
            }
            if (c.color_wig_tracks_by) {
                const col = c.color_wig_tracks_by;
                const vc = this.dataStore.getValueToColor(col);
                //color any wig tracks
                for (const v in vc) {
                    const tr = this.browser.tracks[`${col}|${v}`];
                    if (tr) {
                        tr.config.color = vc[v];
                    }
                }
            }
            if (!this.bamscatrack) {
                const g = this.config.genome_location;
                if (g) {
                    this.browser.update(
                        g.chr,
                        Math.round(g.start),
                        Math.round(g.end),
                        true,
                    );
                } else {
                    this.onDataHighlighted({ indexes: [0] });
                }
            }
            //setup links to other datastores if specified
            //i.e. the browser update if a region is highlighted in another datastore
            this.links=[];
            if (c.sync_with_datastores){
                const dsources = useDataSources();  
                for (let ds of c.sync_with_datastores){
                    const dobj= dsources.find(x=>x.name===ds);
                    if (dobj){
                        this.linkToOtherDataStore(dobj.dataStore);
                    }
                    
                }

            }
        });
    }


    //link this genome browser to another datastore
    async linkToOtherDataStore(dataStore){
        const gb = dataStore.genome_browser;
        if (gb){
            const tname = `${dataStore.name}_base_track`;
            //ensure  the genomic coordinates are loaded in the datastore
            //not sure why loadColumn only loads a single column,
            //but I don't want to mess with it
            for (const col of gb.location_fields){
                await loadColumn(dataStore.name,col);
            }
            //display the the default track from the other datastore
            //in the browser
            if (!this.browser.getTrack(tname)){
                this.browser.addTrack({
                    short_label:gb.default_track.label,
                    track_id:tname,
                    url: gb.default_track.url,
                    decode_function:"generic",
                    displayMode:"EXPANDED",
                    height:20
                },2);
                this.browser.update();
            }
            //listen to the other datastore for highlighted data
            const lid = `gb_${this.config.id}_${dataStore.name}`;
            dataStore.addListener(lid,(type,data)=>{
                if (type==="data_highlighted"){
                    this.onDataHighlighted(data);
                }
            });
            //store info about the link
            this.links.push({
                listener_id:lid,
                dataStore:dataStore
            });
        }
    }

    onBrowserAction(type, data, e) {
        switch (type) {
            case "featureclick":
                if (data.track.config.track_id === "_base_track") {
                    const fIndex = data.feature.data[0];
                    this.dataStore.dataHighlighted(
                        [Number.parseInt(fIndex)],
                        this,
                    );
                }
                if (data.track.config.track_id === "_bam") {
                    console.log(data.feature);
                }
                break;
            case "range_selected":
                if (!e.ctrlKey) {
                    this.browser.update(
                        data["chr"],
                        data["start"],
                        data["end"],
                    );
                } else if (this.bamscatrack) {
                    const ids = this.bamscatrack.getIdsInRange(data);
                    this.cellDim.filter("filterOnIndex", [], ids);
                    this.browser.setHighlightedRegion(data, "_filter", "blue");
                    this.resetButton.style.display = "inline";
                    this.browser.update();
                }
                break;
            case "view_changed":
                this.locText.value = `${data.chr}:${data.start}-${data.end}`;
        }
    }

    //called when chart is created passing the linked datastore
    //and the index key to id as well as the function to get data fron
    //the datastore (ensure column and index are loaded before calling
    //the function)
    setupLinks(dataStore, index, getDataFunction) {
        if (!this.browser) {
            //could change this so it could `await browser`
            //but I don't want to change anything substantial here when
            //I don't have a working example of this chart in any projects.
            setTimeout(
                () => this.setupLinks(dataStore, index, getDataFunction),
                100,
            );
        }

        this.dataLink = { dataStore, index, getDataFunction };
        if (!this.bamscatrack) {
            return;
        }
        const cat = this.config.cluster_reads;

        dataStore.addListener(`gb_${this.config.id}`, (type, data) => {
            if (type === "filtered") {
                this.filterReads(data);
            }
            //has the data used to cluster reads changed?
            //if so update the atac track
            else if (type === "data_changed") {
                const cl = this.config.cluster_reads;
                if (data.columns.indexOf(cl) !== -1) {
                    this.changeClusters(cl, true);
                }
            }
        });
        getDataFunction([cat], () => {
            const ind = dataStore.getColumnIndex(index);
            const colors = dataStore.getColumnColors(cat);
            const col = dataStore.columnIndex[cat];
            this.bamscatrack.setCategories(cat, col.data, col.values, colors);
            //get basic dimension from the other datastore
            this.cellDim = dataStore.getDimension("base_dimension");

            this.bamscatrack.addIndex(ind);

            if (this.config.color_by) {
                this.colorByColumn(this.config.color_by);
            }

            //this.onDataHighlighted({indexes:[0]});
            /*const cf = dataStore.getColorFunction(cat);
            this.browser.setTrackColorFunction("_atac_bam_track",(feature)=>{
                const bc = ("test1#"+feature.tagBA.CB)+"";
                const ci = ind[bc];
                return cf(ci);
            })*/
            const g = this.config.genome_location;
            if (g) {
                this.browser.update(g.chr, g.start, g.end, true);
            } else {
                this.onDataHighlighted({ indexes: [0] });
            }
        });
    }

    changeClusters(column, recalculateCats = false) {
        this.config.cluster_reads = column;
        const ds = this.dataLink.dataStore;
        this.dataLink.getDataFunction([column], () => {
            const colors = ds.getColumnColors(column);
            const col = ds.columnIndex[column];
            this.bamscatrack.setCategories(
                column,
                col.data,
                col.values,
                colors,
                true,
                recalculateCats,
            );
            this.browser.update();
        });
    }
    themeChanged() {
        this.browser.repaint(true, true);
    }

    showHideDialog() {
        const b = this.browser;
        const controls = b.track_order.map((x) => {
            const c = b.tracks[x].config;
            return {
                type: "checkbox",
                id: c.track_id,
                label: c.short_label,
                current_value: !c.hide,
            };
        });
        new CustomDialog({
            title: "Show/Hide Tracks",
            controls: controls,
            doc: this.__doc__,
            buttons: [
                {
                    text: "Update",
                    method: (vals) => {
                        for (const id in vals) {
                            b.tracks[id].config.hide = !vals[id];
                        }
                        b.update();
                    },
                },
            ],
        });
    }

    getConfig() {
        //this method is now called before the browser is created
        //the values it adds to the config should be added dynamically
        //when the user adds a track etc
        const config = super.getConfig();
        const b = this.browser;
        if (b){
            config.genome_location = {
                chr: b.chr,
                start: b.start,
                end: b.end,
            };
            config.tracks = this.browser.getAllTrackConfigs();
        }
    
        return config;
    }

    filterReads(data) {
        if (data === this.cellDim) {
            return;
        }
        if (data === "all_removed") {
            this.bamscatrack.filterReads(null, null);
            this.resetButton.style.display = "none";
            this.browser.removeHighlightedRegion("_filter");
        } else {
            this.bamscatrack.filterReads(
                this.dataLink.dataStore.filterArray,
                this.cellDim.filterArray,
            );
        }
        this.browser.update();
    }
    //called if wig tracks associated with a column
    //obsolete???
    createColumnLinks(dataStore, columns, func) {
        this.dataLink = {
            dataStore: dataStore,
            columns: columns,
            getDataFunction: func,
        };
        //only interested in the first column
        const col = columns[0].col;
        const vc = this.dataStore.getValueToColor(col);
        //color the wig tracks properly
        for (const v in vc) {
            const tr = this.browser.tracks[`${col}|${v}`];
            if (tr) {
                tr.config.color = vc[v];
            }
        }
        // this.onDataHighlighted({indexes:[0]});
    }

    _setPhasedSNPs(info, color) {
        if (info) {
            const d = info.split(",").map((x) => x.split(":"));
            for (const i of d) {
                const st = Number.parseInt(i[1]);
                this.browser.setHighlightedRegion(
                    { chr: i[0], start: st, end: st + 1 },
                    i[0] + i[1],
                    color,
                );
            }
        }
    }

    //called when a feature is highlighted
    onDataHighlighted(data) {
        //if (data.source === this) {
        //    return;
        //}
        const ds = data.dataStore || this.dataStore;
        const p = ds.genome_browser.location_fields;
        let  vm = this.dataStore===ds?this.config.view_margins:ds.genome_browser?.default_parameters?.view_margins;
        vm = vm || this.config.view_margins
        const o = ds.getRowAsObject(data.indexes[0], p);
        //some basic checks
        const st = o[p[2]] > o[p[1]] ? o[p[1]] : o[p[2]];
        let en = o[p[2]] > o[p[1]] ? o[p[2]] : o[p[1]];
        const rowData = data.data;

        //this should be done elsewhere
        const phased = [];
        const not_phased = [];
        if (rowData?.phase_data) {
            this.browser.removeAllHighlightedRegions();
            this._setPhasedSNPs(rowData.more_peak_haplotype, "blue");
            this._setPhasedSNPs(rowData.less_peak_haplotype, "red");
            for (const id of this.browser.track_order) {
                const pd = rowData.phase_data[id];
                if (!id.startsWith("Don")) {
                    continue;
                }
                if (pd) {
                    if (pd.ratio[0] === 0 && pd.ratio[1] === 0) {
                        this.browser.setTrackAttribute(id, "phase_data", null);
                    } else {
                        this.browser.setTrackAttribute(id, "phase_data", {
                            ratio: pd.ratio[0] / (pd.ratio[0] + pd.ratio[1]),
                            st,
                            en,
                        });
                    }
                    phased.push(id);
                } else {
                    this.browser.setTrackAttribute(id, "phase_data", null);
                    not_phased.push(id);
                }
            }
            const new_order = phased.concat(not_phased);
            const track_order = [];
            let index = 0;
            for (const id of this.browser.track_order) {
                if (id.startsWith("Don")) {
                    track_order.push(new_order[index]);
                    index++;
                } else {
                    track_order.push(id);
                }
            }
            //this.browser.track_order=track_order;
        }
        const fcp = this.config.feature_present_column;
        if (fcp) {
            const ids = this.dataStore
                .getRowText(data.indexes[0], fcp)
                .split(", ");
            for (const id of ids) {
                this.browser.setTrackAttribute(id, "color", "#35de26");
            }
            const not = this.dataStore
                .getColumnValues(fcp)
                .slice(0)
                .filter((x) => ids.indexOf(x) === -1);
            for (const id of not) {
                this.browser.setTrackAttribute(id, "color", "#939c92");
            }
            const new_order = ids.concat(not);
            const track_order = [];
            let index = 0;
            for (const id of this.browser.track_order) {
                if (new_order.indexOf(id) !== -1) {
                    track_order.push(new_order[index]);
                    index++;
                } else {
                    track_order.push(id);
                }
            }
            //this.browser.track_order=track_order;
        }

        const textColor = getComputedStyle(this.contentDiv).getPropertyValue('color');
        // set the highlighted region
        this.browser.removeAllHighlightedRegions();
        if (this.config.highlight_selected_region){
            this.browser.setHighlightedRegion({chr:o[p[0]],start:st-1,end:en},"_highlight",textColor,0.3);
        }
        let margin = 1000;
        if (vm.type === "percentage") {
            en = st === en ? st + 1 : en;
            margin = Math.round(en - st) * (vm.value / 100);
        } else {
            margin = vm.value;
        }
        if (vm.type === "fixed_length") {
            const mid = Math.round((en - st) / 2) + st;
            margin = Math.round(margin / 2);
            this.browser.update(o[p[0]], mid - margin, mid + margin);
        } else {
            this.browser.update(o[p[0]], st - margin, en + margin);
        }
    }

    getImage(callback, type) {
        if (type === "svg") {
            this.browser.repaint(true, false, callback);
        } else {
            callback(this.browser.canvas);
        }
    }

    configure_tracks() {
        this.dataSrore;
    }

    getColorOptions() {
        return {
            colorby: "all",
            has_default_color: true,
        };
    }

    _calculatePosition(text) {
        text = text.replace(/,/g, "");

        const arr = text.split(":");
        let chr = null;
        let pos = null;
        if (arr.length === 1) {
            chr = this.browser.getPosition().chr;
            pos = arr[0];
        } else {
            chr = arr[0];
            pos = arr[1];
        }
        const arr2 = pos.split("-");
        return {
            chr: chr,
            start: Number.parseInt(arr2[0]),
            end: Number.parseInt(arr2[1]),
        };
    }

    colorByColumn(column, update = true) {
        const colorFunc = this.getColorFunction(column);
        this.browser.setTrackColorFunction("_base_track", (feature) => {
            return colorFunc(feature.data[0]);
        });
        if (update) {
            this.browser.update();
        }
    }

    setLabelFunction(column) {
        if (!column) {
            this.browser.setTrackLabelFunction("_base_track", null);
            this.config.feature_label = undefined;
        } else {
            this.config.feature_label = column;
            this.browser.setTrackLabelFunction("_base_track", (feature) =>
                this.dataStore.getRowText(feature.data[0], column),
            );
        }
    }

    onDataFiltered(data) {
        if (!this.browser) {
            setTimeout(() => this.onDataFiltered(data), 100);
            return;
        }
        if (!this.dataStore.isFiltered()) {
            this.browser.setTrackFeatureFilter("_base_track", null);
        } else {
            this.browser.setTrackFeatureFilter("_base_track", (feature) => {
                if (this.dataStore.filterArray[feature.data[0]] === 0) {
                    return true;
                }
                return false;
            });
        }

        this.browser.update();
    }

    removeFilter() {
        this.cellDim.removeFilter();
        this.browser.removeHighlightedRegion("_filter");
        this.browser.update();
    }

    changeBaseDocument(doc) {
        this.browser.closeAllDialogs();
        super.changeBaseDocument(doc);
        this.browser.changeBaseDocument(doc);
    }

    setSize(x, y) {
        super.setSize(x, y);
        this.browser.setSize();
        this.browser.update();
    }

    remove(notify = true) {
        if (this.cellDim) {
            this.cellDim.destroy(notify);
            this.dataLink.dataStore.removeListener(`gb_${this.config.id}`);
        }
        //remove any listeners to other datastores
        for (let l of this.links){
            l.dataStore.removeListener(l.listener_id);
        }
        super.remove();
    }

    getSettings() {
        const settings = super.getSettings();
        const cols = this.dataStore.getColumnList();
        const c = this.config;
        cols.push({ name: "None", field: "_none" });
        if (this.bamscatrack) {
            const cats = this.dataLink.dataStore.getColumnList("text");
            settings.push({
                label: "Group By",
                type: "dropdown",
                values: [cats, "name", "field"],
                current_value: c.cluster_reads,
                func: (x) => {
                    this.changeClusters(x);
                },
            });
        }
        settings.push({
            label: "Feature Label",
            type: "dropdown",
            values: [cols, "name", "field"],
            current_value: c.feature_label || "_none",
            func: (x) => {
                x = x === "_none" ? null : x;
                this.setLabelFunction(x);
                this.browser.update();
            },
        });
        settings.push({
            label: "View margin type",
            type: "radiobuttons",
            choices: [
                ["% of feature length", "percentage"],
                ["absolute (bp)", "absolute"],
                ["Fixed Length (bp)", "fixed_length"],
            ],
            current_value: c.view_margins.type,
            func: (x) => {
                c.view_margins.type = x;
                const d = this.dataStore.getHighlightedData();
                if (d) {
                    this.onDataHighlighted({ indexes: [d] });
                }
            },
        });
        settings.push({
            label: "View margin length",
            type: "text",
            current_value: c.view_margins.value,
            only_update_on_enter: true,
            func: (x) => {
                x = Number.parseInt(x);
                x = Number.isNaN(x)
                    ? c.view_margins.type === "percentage"
                        ? 20
                        : 1000
                    : x;
                c.view_margins.value = x;
                const d = this.dataStore.getHighlightedData();
                if (d) {
                    this.onDataHighlighted({ indexes: [d] });
                }
            },
        });
        
        settings.push({
            label:"Highlight Selected Region",
            type:"check",
            current_value:c.highlight_selected_region,
            func:x=>{
                c.highlight_selected_region=x;
                const d = this.dataStore.getHighlightedData();
                if (d){
                    this.onDataHighlighted({indexes:[d]});
                }
              
            }
        })
        return settings;
    }
}

BaseChart.types["genome_browser"] = {
    class: GenomeBrowser,
    name: "Genome Browser",
    methodsUsingColumns: ["setLabelFunction"],
    configEntriesUsingColumns: [
        "feature_label",
        "color_wig_tracks_by",
        "feature_present_column",
    ],

    params: [],
    required: ["genome_browser"],
    init: (config, dataSource) => {
        //set the chr,start,finish columns
        const gb = dataSource.genome_browser;
        config.param = gb.location_fields.slice(0);
        //legacy
        if (gb.feature_label) {
            config.feature_label = gb.feature_label;
        }
        if (gb.default_parameters) {
            Object.assign(config, gb.default_parameters);
        }

        //set the base track corresponding to the datastore
        const df = {
            short_label: gb.default_track.label,
            url: gb.default_track.url,
            track_id: "_base_track",
            decode_function: "generic",
            displayMode: "EXPANDED",
        };
        if (gb.default_track.type) {
            df.type = gb.default_track.type;
        }
        config.tracks = [df];
        if (gb.default_track_parameters) {
            Object.assign(df, gb.default_track_parameters);
        }
        //add the atac track if specified
        if (gb.atac_bam_track) {
            //field to cluster reads
            config.cluster_reads = gb.atac_bam_track.cluster_reads;
            config.tracks.push({
                short_label: "Coverage",
                height: 400,
                track_id: "_atac_bam_track",
                url: gb.atac_bam_track.url,
                type: "bam_sca_track",
            });
        }
        //add any default tracks
        config.default_track = gb.default_track.url;
        if (gb.default_tracks) {
            for (const t of gb.default_tracks) {
                config.tracks.push(JSON.parse(JSON.stringify(t)));
            }
        }
    },
};

export default GenomeBrowser;
