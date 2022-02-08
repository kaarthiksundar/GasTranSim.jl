from attr import validate
from jsonschema import validate
import argparse
from os import listdir
from os.path import isfile, join
import json

def validate_json(data_file, schema_file):
    f = open(data_file)
    data = json.load(f)
    f.close()
    f = open(schema_file)  
    schema = json.load(f)
    f.close()
    validate(instance=data, schema=schema)

parser = argparse.ArgumentParser()
parser.add_argument('--datafolder', type=str, required=True)
args = parser.parse_args()


onlyfiles = [f for f in listdir(args.datafolder) if isfile(join(args.datafolder, f))]

for file in onlyfiles:
    full_filename = args.datafolder + file
    if file.startswith('network'):
        validate_json(full_filename, './case-name/network-schema.json')
    if file.startswith('bc'):
        validate_json(full_filename, './case-name/bc-schema.json')
    if file.startswith('ic'):
        validate_json(full_filename, './case-name/ic-schema.json')
    if file.startswith('params'):
        validate_json(full_filename, './case-name/params-schema.json')
