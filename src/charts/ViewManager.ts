import { action, observable } from "mobx";

class ViewManager {
    @observable accessor current_view;
    @observable accessor all_views;

    constructor(current_view = "", all_views: string[] = []) {
        this.current_view = current_view;
        this.all_views = all_views;
    }

    @action
    setView(view: string) {
        this.current_view = view;
    }

    @action
    setAllViews(all_views: string[]) {
        this.all_views = all_views;
    }

};

export default ViewManager;