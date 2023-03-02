"""
Script to extract biobanks data from rd connect finder source.
It gets metadata from the biobanks and clean some data
"""

import argparse
import json
import os
import re


def clean_orphacodes(code):
    code = code.strip().replace('ORPHA ', 'ORPHA')
    if code in ('sdfsdfs', 'ORPHA', '-'):  # handles dirty data in input
        return []
    if code == 'Beta thalassaemia':  # error in data in input
        return ['ORPHA:848']
    if code == 'FGFR2':  # Error in data in input
        return ['ORPHA:87']
    m = re.search("([A-Z]{2}: ORPHA)[0-9]+", code)  # handles dirty input with a two letters prefix
    if m:
        code = code.replace(m.groups()[0], 'ORPHA')
    if re.match('ORPHA[0-9]+( {2})ORPHA[0-9]+', code):  # handles dirty input with orphacodes separated by two spaces
        code = code.replace('  ', ',')
    if ',' in code:
        return [clean_orphacodes(c)[0] for c in code.split(',') if c != '']
    if ';' in code:
        return [clean_orphacodes(c)[0] for c in code.split(';') if c != '']
    if 'ORPHA' in code:
        return [code.replace('ORPHA', 'ORPHA:')]
    if code not in ('', '*'):
        return [f'ORPHA:{code}']

    return []


def clean_icd10(code):
    code = code.replace('*', '').replace(' ', '').strip()
    # Some strange cases are handled singularly
    if code == 'sdf':
        return []
    if re.match('([A-Z]\.)[0-9]{2}\.[0-9]', code):
        code = code.replace('.', '', 1)
    m = re.match('.*([A-Z][0-9]{2},[0-9])', code)
    if m:
        # handles dirty codes with commas instead of points (Q78,4)
        code = code.replace(m.groups()[0], m.groups()[0].replace(',', '.'))
    if code == "G30-G32":
        return ["urn:miriam:icd:G30-G32"]
    if code == "E71.310-11-12":
        return ['urn:miriam:icd10cm:E71.310', 'urn:miriam:icd10cm:E71.311', 'urn:miriam:icd10cm:E71.312']
    if code in ('G12.0 G12.1 G12.1', 'G12.0 G12.1', 'G12.0G12.1G12.1'):
        return ['urn:miriam:icd:G12.0', 'urn:miriam:icd:G12.1']
    if ',' in code:
        return [clean_icd10(c)[0] for c in code.split(',')]
    if ';' in code:
        return [clean_icd10(c)[0] for c in code.split(';')]
    if code != '':
        if len(code) > 5:
            # it's icd 10 cm
            return [f'urn:miriam:icd10cm:{code}']
        else:
            return [f'urn:miriam:icd:{code}']
    return []


def clean_omim(code):
    code = code.replace('\n', ',').strip()
    if ',' in code:
        return [clean_omim(c)[0] for c in code.split(',') if c != '']
    elif code.strip() != '':
        return [f'OMIM:{code.strip()}']
    return []


def clean_gene(code):
    code = code.strip()
    if code == 'dfsdf':
        return []
    elif ',' in code:
        return [clean_gene(c)[0] for c in code.split(',')]
    elif code != '':
        return [code]
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
