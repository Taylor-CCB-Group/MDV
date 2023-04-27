import re
import json
import pandas as pd
import argparse
import os
import gzip
import shutil
# import numpy as np

parser = argparse.ArgumentParser(
    description='Process a csv table into MDV static site format'
)
parser.add_argument('-i', '--input', help='csv file to process')
parser.add_argument('-o', '--outdir', help='output folder')
parser.add_argument('--discard_redundant', help='discard redundant columns', default=False)

args = parser.parse_args()
filename = args.input
while not os.path.exists(filename):
    input(f'file "{filename}" does not exist. Press enter to try again.')
    filename = input('CSV input file: ')
basename = os.path.basename(filename)
outdir = args.outdir
if not outdir:
    outdir = input('output folder: ')
os.umask(0)
if not os.path.exists(outdir):
    os.makedirs(outdir)

### todo: non-hacky image handling
indir = os.path.dirname(filename)
if os.path.exists(os.path.join(indir, 'images')) and not os.path.exists(os.path.join(outdir, 'images')):
    try:
        shutil.copytree(os.path.join(indir, 'images'), os.path.join(outdir, 'images'))
    except:
        pass

print('reading csv...')
df = pd.read_csv(filename)

types = {
    'float64': 'double',
    'int64': 'integer',
    'O': 'text',
    'object': 'text',
    'bool': 'boolean', # 0 or 1? currently breaks this script if I make this integer
}

def get_column_type(name):
    v = df[name]
    unique_values = len(set(v))
    type = types[str(v.dtype)]
    if type == 'text' and unique_values == v.size:
        return 'unique'
    if type is None:
        raise ValueError(f'unknown type {v.dtype} for {name}')
    return type


def get_quantiles(col):
    qs = {}
    for q in ["0.001", "0.01", "0.05"]:
        q1 = col.quantile(float(q))
        q2 = col.quantile(1-float(q))
        qs[q] = [q1, q2]
    return qs

def get_text_indices(col):
    values = list(set(col))
    val_dict = {value: i for i, value in enumerate(values)}
    return [val_dict[v] for v in col], values

def get_column_groups():
    cols = {}
    for name in df.columns:
        m = re.search(r'(.+)_(\d+)(.*)', name)
        if not m:
            continue
        group_name = f'{m.group(1)}_{m.group(3)}'
        if group_name not in cols:
            cols[group_name] = {'name': group_name, 'columns': []}
        # num = int(m.group(2))
        cols[group_name]['columns'].append(name)
    return [cols[k] for k in cols]

def get_datasource():
    '''
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
    '''
    descriptor = {
        "name": basename,
        "size": df.shape[0],
        "images": {
            "images": {
                "base_url": "./images/",
                "type": "png",
                "key_column": "Index"
            }
        },
        "columns": []
    }
    for name in df.columns:
        col = df[name]
        if args.discard_redundant and len(set(col)) == 1:
            df.drop(name, axis=1, inplace=True)
            continue
        datatype = get_column_type(name)
        col_desc = { "datatype": datatype, "name": name, "field": name }
        if datatype == 'double' or datatype == 'integer':
            col_desc['minMax'] = [min(col), max(col)]
            col_desc['quantiles'] = get_quantiles(col)
        elif datatype == 'text':
            indices, values = get_text_indices(col)
            col_desc['values'] = values
            df[name] = indices
        elif datatype == 'unique':
            col_desc['stringLength'] = max([len(v) for v in col])
        descriptor['columns'].append(col_desc)
    descriptor['columnGroups'] = get_column_groups()
    return descriptor

def replace_text_values(col, values):
    val_dict = {value: i for i, value in enumerate(values)}
    return [val_dict[v] for v in col]

def get_views():
    return {basename: {"name": basename, 'initialCharts': {basename: []}}}

def get_state():
    return {"all_views": [basename], "initial_view": basename}

def convert_data_to_binary(df):
    '''
    Converts the dataframe to binary format.
    '''
    dfile = f'{outdir}/{basename}.b'
    o = open(dfile, 'wb')
    index = {}
    current_pos = 0
    for name in df.columns:
        # 'integer' and 'double' should be converted to float32 according to the spec
        type = get_column_type(name)
        if type == 'integer' or type == 'double':
            df[name] = df[name].astype('float32')
        comp = gzip.compress(df[name].to_numpy().tobytes())
        new_pos = current_pos + len(comp)
        index[name] = [current_pos, new_pos-1]
        o.write(comp)
        current_pos = new_pos
    o.close()
    ifile = dfile[:dfile.rindex('.')] + '.json'
    with open(ifile, 'w') as f:
        f.write(json.dumps(index))

def main():
    if not os.path.exists(outdir):
        print('creating output directory')
        os.makedirs(outdir)

    with open(f'{outdir}/datasources.json', 'w') as f:
        print('writing datasources.json')
        f.write(json.dumps([get_datasource()]))

    print('writing data binary')
    convert_data_to_binary(df)

    with open(f'{outdir}/views.json', 'w') as f:
        f.write(json.dumps(get_views()))
    with open(f'{outdir}/state.json', 'w') as f:
        f.write(json.dumps(get_state()))

if __name__ == '__main__':
    main()