import { getRandomString } from "@/utilities/Utilities";
import type BaseChart from "./BaseChart";
import { makeAutoObservable } from "mobx";


export function serialiseConfig<T extends BaseChart<any>>(chart: T) {
    const { config } = chart;
    // thinking about config vs state, in the context of dynamic virtual columns...
    // if we just have a record of paramSpecs, then we can use those in this representation
    // while the config object itself will have the actual evaluated runtime state of the values
    // then if we also keep a record of other `configEntriesUsingColumns`, that should, hypothetically,
    // be enough information?
    // - wouldn't like to count on that.
    // If we make sure that any special config values are of some type other than string, 
    // then we should be easily able to traverse the config object for them... however, it works the other way around:
    // the config object is liable to have computed strings & we want to get the 'special value' from that.

    // get the BaseChart.types entry for this chart, 
    // use it to determine any `configEntriesUsingColumns`...
    const serialized = JSON.parse(JSON.stringify(config));
    console.log('processed config:', serialized);
    return serialized;
}

export function initialiseConfig<T extends BaseChart<any>>(originalConfig: any, chart: T) {
    let config = JSON.parse(JSON.stringify(originalConfig));
    if (!config.id) {
        // what about when we duplicate a chart?
        config.id = getRandomString();
    }

    // makeAutoObservable(config);
    Object.defineProperty(chart, "config", {
        get: () => config,
        set: (v) => {
            // this will be invoked after initialiseConfig returns in BaseChart constructure,
            // so we no long makeAutoObservable(config) above, but this will apply it there
            config = v; // re-assigning config ref in this closure is not really relevant now;
            //we are using makeAutoObservable and not referring to config directly...
            makeAutoObservable(config);
        },
    });

    setTimeout(() => {
        // defer this until after the constructor has finished
        chart.mobxAutorun(() => {
            // setTitle also sets config.title - which is what we're observing here...
            // so we got warnings about setTitle mutating config.title 'outside an action' (at least it didn't go into infinite loop).
            // perhaps title.textContent could be a computed value and we may not need this autorun at all...
            // this.title.textContent = config.title;
            // we can now safely call this.setTitle() without warnings as it avoids unnecessary config.title changes.
            chart.setTitle(config.title);
        });
    }, 0);
    // note: a previous version of this used makeObservable for keeping track of onDataFiltered...
    // that worked, with extra extraneous number that changed to be observed by the hook...
    // What I have now done is change DataStore to be observable, and added a method for getting filtered indices
    // in a way that can be shared by different charts (react or otherwise).


    return config;
}

