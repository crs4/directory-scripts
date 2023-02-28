"""
Script to extract biobanks data from rd connect finder source.
It gets metadata from the biobanks and clean some data
"""

import argparse
import json
import os


def clean_orphacodes(code):
    if ',' in code:
        return [clean_orphacodes(c)[0] for c in code.split(',')]
    elif 'ORPHA' in code:
        return [code.strip().replace('ORPHA', 'ORPHA:')]
    elif code not in ('', '*'):
        return [f'ORPHA:{code.strip()}']
    else:
        return []


def clean_icd10(code):
    if ',' in code:
        return [clean_icd10(c)[0] for c in code.split(',')]
    elif ';' in code:
        return [clean_icd10(c)[0] for c in code.split(';')]
    elif code.strip() != '':
        if len(code) > 5:
            # it's icd 10 cm
            return [f'urn:miriam:icd10cm:{code.strip()}']
        else:
            return [f'urn:miriam:icd:{code.strip()}']
    else:
        return []


def clean_omim(code):
    code = code.replace('\n', ',')
    if ',' in code:
        print(code.strip())
        return[clean_omim(c)[0] for c in code.split(',') if c != '']
    elif code.strip() != '':
        return [f'OMIM:{code.strip()}']
    return []


def clean_gene(code):
    if ',' in code:
        return [clean_gene(c)[0] for c in code.split(',')]
    elif code.strip() != '':
        return [code.strip()]
    return []


def clean_material(code):
    if ',' in code:
        return [clean_material(c.strip())[0] for c in code.split(',')]
    elif ';' in code:
        return [clean_material(c.strip())[0] for c in code.split(';')]
    elif ' - ' in code:
        return [clean_material(c.strip())[0] for c in code.split('-')]
    elif code.strip() == '':
        return []
    return [code.strip()]


def file_exist(file_argument):
    if os.path.exists(file_argument):
        return file_argument
    raise argparse.ArgumentTypeError("File {} does not exist".format(file_argument))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file', '-i', dest='rd_connect_file', type=file_exist, required=True,
                        help='The JSON input file with data of RD Connect\'s biobanks')
    args = parser.parse_args()

    with open(args.rd_connect_file) as f:
        rd_connect_data = json.load(f)

    rdc_biobanks = [b for b in rd_connect_data['allData'] if b['type'] in 'biobank']

    for b in rdc_biobanks:
        for d in b['diseases']:
            d['orphacode'] = clean_orphacodes(d['orphacode'])
            d['icd10'] = clean_icd10(d['icd10'])
            d['omim'] = clean_omim(d['omim'])
            d['gene'] = clean_gene(d['gene'])
        b['bb_core']['Additional_Biomaterial_available'] = clean_material(
            b['bb_core']['Additional_Biomaterial_available'])
        b['bb_core']['Biomaterial_Available'] = json.loads(b['bb_core']['Biomaterial_Available'])

    with open('rd_connect_biobanks.json', 'w') as o:
        json.dump({'allData': rdc_biobanks}, o, indent=2)
