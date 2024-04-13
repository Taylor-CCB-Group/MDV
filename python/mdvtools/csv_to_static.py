import re
import json
import pandas as pd
import argparse
import os
import gzip
import shutil
# import numpy as np

parser = argparse.ArgumentParser(
    description="Process a csv table into MDV static site format"
)
parser.add_argument("-i", "--input", help="csv file to process", default="data.csv")
parser.add_argument("-o", "--outdir", help="output folder")
parser.add_argument(
    "--discard_redundant", help="discard redundant columns", default=False
)
# grouping options? rename columns so numbers have leading zeros?
parser.add_argument("-g", "--group", help="group columns by regex", default=False)
parser.add_argument(
    "--group_by", help="group columns by regex", default=r"(.*?)(\d+)(.*)"
)
parser.add_argument("-s", "--separator", help="multitext separator", default=";")
parse_multitext = True

args = parser.parse_args()
separator = args.separator
filename = args.input
while not os.path.exists(filename):
    input(f'file "{filename}" does not exist. Press enter to try again.')
    filename = input("CSV input file: ")
basename = os.path.basename(filename)
outdir = args.outdir
if not outdir:
    outdir = input("output folder: ")
os.umask(0)
if not os.path.exists(outdir):
    os.makedirs(outdir)

### todo: non-hacky image handling
indir = os.path.dirname(filename)
has_images = os.path.exists(os.path.join(indir, "images"))
if has_images and not os.path.exists(os.path.join(outdir, "images")):
    try:
        shutil.copytree(os.path.join(indir, "images"), os.path.join(outdir, "images"))
    except Exception:
        pass

print("reading csv...")
df = pd.read_csv(filename)

types = {
    "float64": "double",
    "int64": "integer",
    "O": "text",
    "object": "text",
    "bool": "boolean",  # 0 or 1? currently breaks this script if I make this integer
}

# if we were processing multiple sources, we should review global df / col_types...
col_types = {}


def rename_columns():
    return
    # rename columns that are numbers to have leading zeros
    # this is so they sort correctly in the UI
    for name in df.columns:
        m = re.search(args.groub_by, name)
        if not m:
            continue
        new_name = f"{m.group(1)}_{m.group(2).zfill(3)}{m.group(3)}"
        df.rename(columns={name: new_name}, inplace=True)


def get_column_type(name):
    # get_column_type is called from get_datasource() then convert_data_to_binary()
    # second call was getting wrong type, so remembering the values should help.
    if name in col_types:
        return col_types[name]
    v = df[name]
    unique_values = set(v)
    dtype = str(v.dtype)
    ttype = types[dtype]
    # if dtype == 'text' and len(unique_values) == v.size:
    #     print(f'unique text column "{name}"')
    #     ttype = 'unique'
    if ttype == "text" and parse_multitext:
        # does it look like comma-separated tags?
        # 'argument of type 'bool' is not iterable'???
        # when we have something like "unique_values: {False, True, nan}"
        # print(f'{name}: ({type}) unique_values: {unique_values}')
        n = len(unique_values)
        if n > 65536:
            print(f'detected unique column "{name}" (not well tested with this script)')
            ttype = "unique"
        elif n > 256 or any([separator in str(s) for s in unique_values]):
            print(f'detected multitext column "{name}"')
            ttype = "multitext"
    if ttype is None:
        raise ValueError(f"unknown type {v.dtype} for {name}")
    col_types[name] = ttype
    return ttype


def get_quantiles(col):
    qs = {}
    for q in ["0.001", "0.01", "0.05"]:
        q1 = col.quantile(float(q))
        q2 = col.quantile(1 - float(q))
        qs[q] = [q1, q2]
    return qs


def get_text_indices(col):
    values = list(set(col))
    val_dict = {value: i for i, value in enumerate(values)}
    return [val_dict[v] for v in col], [str(s) for s in values]


def get_column_groups():
    col_groups = {}
    for name in df.columns:
        m = re.search(args.group_by, name)
        if not m:
            continue
        group_name = f"{m.group(1)}_{m.group(3)}"
        if group_name not in col_groups:
            col_groups[group_name] = {"name": group_name, "columns": []}
        # num = int(m.group(2))
        col_groups[group_name]["columns"].append(name)
    return [col_groups[k] for k in col_groups]


def get_datasource():
    """
    Has some side effects on the dataframe:
    if args.discard_redundant:
        - removes columns that are redundant (all the same value)
    text columns are converted to indices.

    Outputs a descriptor like this:
    {
        "name": "metric_table",
        "size": number of rows,
        "images": {
            "images": {
                "base_url": "./images/",
                "type": "png",
                "key_column": "image_id"
            }
        }
        "columns": [
            {
                "datatype": "float" | "integer" | "text" | "unique",
                "name": "column_name",
                "field": "column_name",
                "minMax"?: [min, max],
                "quantiles"?: ...,
                "values"?: ['a', 'b', 'c'],
            }
        ]
    }
    """
    descriptor = {"name": basename, "size": df.shape[0], "columns": []}
    if has_images:
        # todo: make this able to take some config, set proper type / key_column etc.
        # ideally, find a column that has values corresponding to the names of images in the folder...
        descriptor["images"] = {
            "images": {"base_url": "./images/", "type": "png", "key_column": "Index"}
        }
    for name in df.columns:
        col = df[name]
        if args.discard_redundant and len(set(col)) == 1:
            df.drop(name, axis=1, inplace=True)
            continue
        datatype = get_column_type(name)
        col_desc = {"datatype": datatype, "name": name, "field": name}
        if datatype == "boolean":
            # would be better to have a separate boolean type
            print(f"converting boolean {name} to number")
            col_desc["datatype"] = "integer"
            col_desc["minMax"] = [0, 1]
        elif datatype == "double" or datatype == "integer":
            col_desc["minMax"] = [min(col), max(col)]
            col_desc["quantiles"] = get_quantiles(col)
        elif datatype == "text" or datatype == "multitext":
            # would be better to have a separate boolean type
            # col_desc['datatype'] = 'text'
            indices, values = get_text_indices(col)
            col_desc["values"] = values
            if datatype == "multitext":
                col_desc["separator"] = separator
            # mutating df here...
            df[name] = indices
        elif datatype == "unique":
            col_desc["stringLength"] = max([len(v) for v in col])
        descriptor["columns"].append(col_desc)
    descriptor["columnGroups"] = get_column_groups()
    return descriptor


def replace_text_values(col, values):
    val_dict = {value: i for i, value in enumerate(values)}
    return [val_dict[v] for v in col]


def get_views():
    return {basename: {"name": basename, "initialCharts": {basename: []}}}


def get_state():
    return {"all_views": [basename], "initial_view": basename}


def convert_data_to_binary(df):
    """
    Converts the dataframe to binary format.
    """
    dfile = f"{outdir}/{basename}.gz"
    o = open(dfile, "wb")
    index = {}
    current_pos = 0
    for name in df.columns:
        # 'integer' and 'double' should be converted to float32 according to the spec
        type = get_column_type(name)
        if type == "integer" or type == "double" or type == "boolean":
            print(f"converting {name} {type} to float32")
            df[name] = df[name].astype("float32")
        if type == "text":
            print(f"converting {name} {type} to uint8")
            df[name] = df[name].astype("uint8")
        if type == "multitext":
            print(f"converting {name} {type} to uint16")
            df[name] = df[name].astype("uint16")
        comp = gzip.compress(df[name].to_numpy().tobytes())
        new_pos = current_pos + len(comp)
        index[name] = [current_pos, new_pos - 1]
        o.write(comp)
        current_pos = new_pos
    o.close()
    ifile = dfile[: dfile.rindex(".")] + ".json"
    with open(ifile, "w") as f:
        f.write(json.dumps(index))


def main():
    rename_columns()
    if not os.path.exists(outdir):
        print("creating output directory")
        os.makedirs(outdir)

    ds = get_datasource()
    with open(f"{outdir}/datasources.json", "w") as f:
        print("writing datasources.json")
        f.write(json.dumps([ds]))

    print("writing data binary")
    convert_data_to_binary(df)

    with open(f"{outdir}/views.json", "w") as f:
        f.write(json.dumps(get_views()))
    with open(f"{outdir}/state.json", "w") as f:
        f.write(json.dumps(get_state()))


if __name__ == "__main__":
    main()
