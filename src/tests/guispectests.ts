import type { GuiSpec, GuiSpecs } from "@/charts/charts";
import { g } from "@/lib/utils";

const untyped = {
    type: 'dropdown',
    label: 'label',
    current_value: 'current_value',
    //@ts-expect-error no implicit any, and we are in an anonymous object
    func: (v) => { },
    values: [[
        { a: '0a', b: '0b' },
        { a: '1a', b: '1b' },
    ], 'a', 'b']
};
// const t = g(untyped);

const a: GuiSpec<'dropdown'> = {
    type: 'dropdown',
    label: 'label',
    //@ts-expect-error
    name: 'name',
    current_value: 'current_value',
    func: (v) => {},
    values: [[
        {a: '0a', b: '0b'},
        {a: '1a', b: '1b'},
    ], 'a', 'b']
}
if (a.values) {
    // if we were really clever we might know that x.a is a string here
    // a.values[0].map(x => x.a);
}

const b = g({
    type: 'dropdown',
    label: 'label',
    current_value: 'current_value',
    //correctly inferred, explicit type is not needed and is be compiler-checked
    func: (v: string) => {},
    values: [['a', 'b']]
})

const c = g({
    type: "check",
    label: "label",
    current_value: true,
    func(v) {
        c.current_value = v;
    },
})
// g hasn't made it know that func is defined - that is correct as it could change
// it has inferred the type correctly, though
c.func?.(false);

const specArray: GuiSpecs = [
    g({
        type: "check",
        label: "label",
        current_value: true,
        func(v) {
            c.current_value = v;
        }
    }),
    g({
        type: "dropdown",
        label: "label",
        current_value: "current_value",
        func(v) {
            b.current_value = v;
            this.current_value = v;
        }
    })
]

specArray.concat([
    g({
        type: "check",
        label: "label",
        current_value: true,
        async func(v) {
            this.current_value = v;
        }
    }),
    g({
        type: "multidropdown",
        label: "label",
        current_value: ["current_value"],
        func(v) {
            this.current_value = v.slice();
        }
    })
]);
