import _ from "lodash";

const useDiff = (input1: object, input2: object) => {
    const diff = (obj1: any, obj2: any) => {
        return _.reduce(
            obj1,
            (result: any, value: any, key: string) => {
                if (_.isPlainObject(value)) {
                    const d = diff(value, obj2[key]);
                    if (!_.isEmpty(d)) {
                        result[key] = d;
                    }
                } else if (!_.isEqual(value, obj2[key])) {
                    // different from SO answer, which wouldn't diff arrays.
                    result[key] = diff(value, obj2[key]);
                } else {
                    delete result[key];
                }
                return result;
            },
            {},
        );
    };

    return diff(input1, input2);
};

export default useDiff;
