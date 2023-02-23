"""
This script gets metadata of Biobanks from RDConnect Finder stored as JSON and
updates the EMX file of the BBMRI Directory creating records for the missing biobanks or updating the ones already
present in the directory.
To get the updated disease codes, the Sample Catalogue is also queried.
"""
import argparse
import json
import math
import os

import pandas as pd
from molgenis.client import Session

BIOBANKS_MAPPING = {
    173631: {
        'biobank': 'bbmri-eric:ID:IT_1384375808568138',
        'collection': 'bbmri-eric:ID:IT_1384375808568138:collection:794f03ad818541b'
    },
    '215190': {
        'biobank': 'bbmri-eric:ID:IT_1539241927435687',
        'collection': 'bbmri-eric:ID:IT_1539241927435687:collection:c9346d1b70d14fc'
    },
    '45274': {
        'biobank': 'bbmri-eric:ID:IT_1383569594559774',
        'collection': 'bbmri-eric:ID:IT_1383569594559774:collection:1444717432314364'
    },
    '77489': {
        'biobank': 'bbmri-eric:ID:IT_1383047508168267',
        'collection': 'bbmri-eric:ID:IT_1383047508168267:collection:1444717360015607'
    }
}

COUNTRIES_CODES = {
    'Germany': 'DE',
    'Hungary': 'HU',
    'Israel': 'IL',
    'Italy': 'IT',
    'Malta': 'MT',
    'Slovenia': 'SI',
    'Spain': 'ES',
    'Turkey': 'TR',
    'United Kingdom': 'UK'
}

MATERIAL_TYPES = {
    'BLOOD': 'OTHER',
    'CELLS': 'OTHER',
    'DNA': 'DNA',
    'LYMPHOBLASTOID CELL LINES': 'CELL_LINES',
    'NON-INMORTALIZED CELL LINES': 'CELL_LINES',
    'OTHER (PLEASE SPECIFY)': 'OTHER',
    'OTHER FLUIDS': 'OTHER',
    'PLASMA': 'PLASMA',
    'PRIMARY CELL LINES': 'CELL_LINES',
    'RNA': 'RNA',
    'SALIVA': 'SALIVA',
    'SERA': 'SERUM',
    'SERUM': 'SERUM',
    'TISSUES': 'TISSUE_PARAFFIN_EMBEDDED',
    'URINE': 'URINE',
    'WHOLE BLOOD': 'WHOLE_BLOOD'
}

RD_NETWORK = 'bbmri-eric:networkID:EU_BBMRI-ERIC:networks:RD-Biobanks'

BIOBANKS_SHEET = 'eu_bbmri_eric_biobanks'
COLLECTIONS_SHEET = 'eu_bbmri_eric_collections'
PERSONS_SHEET = 'eu_bbmri_eric_persons'
ALSO_KNOWN_SHEET = 'eu_bbmri_eric_also_known_in'


class RDConnectImporter:

    def __init__(self, bb_dir_data, rdc_finder_data, sc_url, sc_user, sc_pass):
        self.eric_data = bb_dir_data
        self.rdc_finder_data = rdc_finder_data
        self.sc_url = sc_url
        self.sc_user = sc_user
        self.sc_pass = sc_pass
        self.sc_session = Session(self.sc_url)
        self.sc_session.login(self.sc_user, self.sc_pass)

    def get_country_code(self, rd_connect_country):
        return COUNTRIES_CODES[rd_connect_country]

    def biobank_exist(self, rd_biobank_id):
        return rd_biobank_id in BIOBANKS_MAPPING

    def get_biobank_id(self, rd_biobank_id, biobank_country):
        try:
            return BIOBANKS_MAPPING[rd_biobank_id]['biobank']
        except KeyError:
            return f'bbmri-eric:ID:RD_{biobank_country}:{rd_biobank_id}'

    def get_collection_id(self, rd_biobank_id, biobank_country):
        try:
            return BIOBANKS_MAPPING[rd_biobank_id]['collection']
        except KeyError:
            return f'{self.get_biobank_id(rd_biobank_id, biobank_country)}:collection:MainCollection'

    def generate_contact_id(self, rd_biobank_id, biobank_country):
        return f'bbmri-eric:contactID:RD_{biobank_country}_{rd_biobank_id}'

    def get_material(self, rd_material):
        try:
            return MATERIAL_TYPES[rd_material.upper()]
        except KeyError:
            return None

    def get_donors_and_diseases(self, rd_data):
        def _get_orphacodes(code):
            if ',' in code:
                return [_get_orphacodes(c.strip())[0] for c in code.split(',')]
            elif 'ORPHA' in code:
                return [code.replace('ORPHA', 'ORPHA:')]
            elif code not in ('', '*'):
                return [f'ORPHA:{code}']
            else:
                return []

        num_donors = 0
        diseases = set()
        for d in rd_data['diseases']:
            diseases.update(_get_orphacodes(d['orphacode']))
            if d['icd10'].strip() != '':
                diseases.add(f"urn:miriam:icd:{d['icd10']}")
            num_donors += int(d['number'])

        diseases.update(self.get_diseases_from_sample_catalogue(rd_data['OrganizationID']))
        return num_donors, diseases

    def get_diseases_from_sample_catalogue(self, rd_biobank_id):
        print("getting diseases from sample catalogue")
        diseases = set()
        sc_res = self.sc_session.get('rd_connect_Sample', q=f'BiobankID=={rd_biobank_id}', attributes='Disease')
        if 'errors' not in sc_res:
            print("Done getting additional diseases")
        for r in sc_res:
            diseases.update(
                [d['ID'].replace('urn:miriam:orphanet:', 'ORPHA:') for d in r['Disease'] if "ncit" not in d['ID']])
        return diseases

    def get_organization_info(self, rd_data):
        country_code = self.get_country_code(rd_data['address']['country'])
        biobank_id = self.get_biobank_id(rd_data['OrganizationID'], country_code)
        df = self.eric_data[BIOBANKS_SHEET]
        if df[df.id == biobank_id].empty:
            new_biobank = pd.DataFrame({
                'id': [biobank_id],
                'pid': [biobank_id],
                'name': [rd_data['name']],
                'acronym': [rd_data['bb_core']['acronym']],
                'description': [rd_data['bb_core']['Description']],
                'url': [rd_data['url'][0]],
                'location': [''],
                'country': [country_code],
                'head': [''],
                'contact': [self.generate_contact_id(biobank_id, country_code)],
                'juridical_person': [rd_data['address']['name of host institution']],
                'network': [RD_NETWORK],
                'also_known': [''],
                'collections': [f'{biobank_id}:collection:MainCollection'],
                'capabilities': [''],
                'quality': [''],
                'collaboration_commercial': [''],
                'collaboration_non_for_profit': [''],
                'national_node': [country_code],
                'withdrawn': [True]
            }, index=None)
            self.eric_data[BIOBANKS_SHEET] = pd.concat([self.eric_data[BIOBANKS_SHEET], new_biobank])
        else:
            df.at[df.id == biobank_id, 'network'] = [f"{df.loc[df.id == biobank_id]['network'].values[0]},{RD_NETWORK}"]

    def _concat_cell(self, current_cell, new_string):
        if str(current_cell.values[0]) == 'nan':
            return new_string
        else:
            return f'{current_cell.values[0]},{new_string}'

    def get_collection_info(self, rd_data):
        country_code = self.get_country_code(rd_data['address']['country'])
        rd_biobank_id = rd_data['OrganizationID']
        collection_id = self.get_collection_id(rd_biobank_id, country_code)

        additional = [am.strip() for am in rd_data['bb_core']['Additional_Biomaterial_available'].split(',') if
                      am != '']
        materials = [self.get_material(m) for m in json.loads(rd_data['bb_core']['Biomaterial_Available']) + additional
                     if self.get_material(m) is not None]

        num_of_donors, diseases = self.get_donors_and_diseases(rd_data)

        df = self.eric_data[COLLECTIONS_SHEET]
        if df[df.id == collection_id].empty:
            new_collection = pd.DataFrame({
                'id': [self.get_collection_id(collection_id, country_code)],
                'name': ['Main Collection'], 'acronym': [''], 'description': [''], 'url': [''], 'location': [''],
                'country': [country_code], 'head': [''],
                'contact': [self.generate_contact_id(rd_biobank_id, country_code)],
                'withdrawn': [''], 'national_node': [''],  # TODO
                'parent_collection': [''], 'sub_collections': [''],
                'biobank': [collection_id], 'biobank_label': [rd_data['name']], 'network': [RD_NETWORK],
                'combined_network': [RD_NETWORK], 'also_known': [f'rdconnect:{rd_biobank_id}'],
                'type': ['RD'],  # statically add only rare disease type
                'data_categories': ['BIOLOGICAL_SAMPLES,OTHER'],
                'order_of_magnitude': ['0'],  # TODO: check if it can be retrieved in the samples catalogu]e
                'size': [''],  # TODO: check if it can be retrieved in the samples catalogue
                'categories': [''], 'timestamp': [''], 'quality': [''], 'combined_quality': [''],
                'number_of_donors': [num_of_donors],
                'order_of_magnitude_donors': [int(math.log10(max(1, num_of_donors)))],
                'sex': [''],  # TODO: check if it can be retrieved in the samples catalogue
                'diagnosis_available': [','.join(sorted(diseases))], 'age_low': [''],
                'age_high': [''], 'age_unit': [''], 'materials': [','.join(sorted(materials))],
                'storage_temperatures': [''], 'body_part_examined': [''], 'imaging_modality': [''],
                'image_dataset_type': [''], 'collaboration_commercial': [''], 'collaboration_non_for_profit': [''],
                'data_use': [''], 'commercial_use': [''], 'access_fee': [''], 'access_joint_project': [''],
                'access_description': [''], 'access_uri': [''], 'sop': ['']
            })
            self.eric_data[COLLECTIONS_SHEET] = pd.concat([df, new_collection])
        else:
            df.at[df.id == collection_id, ['also_known', 'network', 'combined_network', 'diagnosis_available']] = \
                [self._concat_cell(df.loc[df.id == collection_id, 'also_known'], f'rdconnect:{rd_biobank_id}'),
                 self._concat_cell(df.loc[df.id == collection_id, 'network'], RD_NETWORK),
                 self._concat_cell(df.loc[df.id == collection_id, 'combined_network'], RD_NETWORK),
                 self._concat_cell(df.loc[df.id == collection_id, 'diagnosis_available'], ','.join(sorted(diseases)))]

    def get_contact_info(self, rd_data):
        country_code = self.get_country_code(rd_data['address']['country'])
        biobank_id = self.get_biobank_id(rd_data['OrganizationID'], country_code)
        return {
            'id': self.generate_contact_id(rd_data['OrganizationID'], country_code),
            'title_before_name': '',
            'first_name': rd_data['main contact']['first name'],
            'last_name': rd_data['main contact']['last name'],
            'title_after_name': '',
            'email': rd_data['main contact']['email'],
            'phone': '',
            'address': '',
            'zip': '',
            'city': '',
            'country': country_code,
            'role': '',
            'biobanks': biobank_id,
            'collections': self.get_collection_id(biobank_id),
            'networks': '',
            'national_node': country_code,
            'withdrawn': True
        }

    def get_also_known_in_info(self, rd_data):
        rd_biobank_id = rd_data['OrganizationID']
        return {
            'id': f'rdconnect:{rd_biobank_id}',
            'name_system': 'RD Connect',
            'pid': rd_biobank_id,
            'url': '',
            'national_node': self.get_country_code(rd_data['address']['country']),
            'withdrawn': ''
        }

    def run(self):
        # eric_data = defaultdict(list)
        for b in self.rdc_finder_data:
            print("Converting biobank: ", b['OrganizationID'])
            print("Getting biobank data: ", b['OrganizationID'])
            self.get_organization_info(b)
            print("Getting collection data: ", b['OrganizationID'])
            self.get_collection_info(b)
            # print("Getting person data: ", b['OrganizationID'])
            # self.get_contact_info(b)
            # print("Getting also known data: ", b['OrganizationID'])
            # self.get_also_known_in_info(b)
            print()
        return self.eric_data


def write_excel(data, output_file):
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        for sheet, df in data.items():
            df.to_excel(writer, sheet_name=sheet, index=False)


def file_exist(file_argument):
    if os.path.exists(file_argument):
        return file_argument
    raise argparse.ArgumentTypeError("File {} does not exist".format(file_argument))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory-emx-file', '-e', dest='directory_file',
                        type=file_exist, required=True,
                        help='Input directory emx file with. A copy with new data will be created')
    parser.add_argument('--rd-connect-input-file', '-r', dest='rd_connect_file', type=file_exist, required=True,
                        help='The JSON input file with data of RD Connect\'s biobanks')
    parser.add_argument('--output', '-o', dest='output_file', type=str, required=False,
                        help='Output EMX file containing updated data', default='bbmri-directory.xlsx')
    parser.add_argument('--sc-url', '-H', dest='sc_url', type=str, required=True,
                        help='Base URL of the RD Connect Samples Catalogue')
    parser.add_argument('--sc-user', '-u', dest='sc_user', type=str, required=True,
                        help='Samples Catalogue user name')
    parser.add_argument('--sc-password', '-p', dest='sc_pwd', type=str, required=True,
                        help='Samples Catalogue user password')

    args = parser.parse_args()

    with open(args.rd_connect_file) as f:
        rd_connect_data = json.load(f)

    rdc_biobanks = [b for b in rd_connect_data['allData'] if
                    b['type'] in 'biobank' and b['OrganizationID'] in (44001, 173631)]
    d_emx = pd.read_excel(args.directory_file, sheet_name=None, engine='openpyxl')

    importer = RDConnectImporter(d_emx, rdc_biobanks, args.sc_url, args.sc_user, args.sc_pwd)
    res = importer.run()
    write_excel(res, args.output_file)
