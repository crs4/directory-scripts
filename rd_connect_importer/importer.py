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
from bs4 import BeautifulSoup
from molgenis.client import Session

ITALIAN_BIOBANKS = [45274, 45401, 62316, 76957, 77088, 77219, 77350, 77489, 77630, 77761, 173631, 215190, 261780]

BIOBANKS_MAPPING = {
    45274: {
        'biobank': 'bbmri-eric:ID:IT_1383569594559774',
        'collection': 'bbmri-eric:ID:IT_1383569594559774:collection:1444717432314364'
    },
    45401: {
        'biobank': 'bbmri-eric:ID:IT_1383929642991480',
        'collection': 'bbmri-eric:ID:IT_1383929642991480:collection:1444717553513748'
    },
    62316: {
        'biobank': 'bbmri-eric:ID:IT_1382965524316631',
        'collection': 'bbmri-eric:ID:IT_1382965524316631:collection:1444717339490516'
    },
    77219: {
        'biobank': 'bbmri-eric:ID:IT_1383230735994332',
        'collection': 'bbmri-eric:ID:IT_1383230735994332:collection:1450446268862465'
    },
    77350: {
        'biobank': 'bbmri-eric:ID:IT_1382965524316631',
        'collection': 'bbmri-eric:ID:IT_1382965524316631:collection:1444717339490516'
    },
    77489: {
        'biobank': 'bbmri-eric:ID:IT_1383047508168267',
        'collection': 'bbmri-eric:ID:IT_1383047508168267:collection:1444717360015607'
    },
    77630: {
        'biobank': 'bbmri-eric:ID:IT_1385652938842205',
        # 'collection': 'bbmri-eric:ID:IT_1385652938842205:collection:1444717851150981'
    },
    173631: {
        'biobank': 'bbmri-eric:ID:IT_1384375808568138',
        'collection': 'bbmri-eric:ID:IT_1384375808568138:collection:794f03ad818541b'
    },
    215190: {
        'biobank': 'bbmri-eric:ID:IT_1539241927435687',
        'collection': 'bbmri-eric:ID:IT_1539241927435687:collection:c9346d1b70d14fc'
    },
    261780: {
        'biobank': 'bbmri-eric:ID:IT_1385652938842205',
        # 'collection': 'bbmri-eric:ID:IT_1385652938842205:collection:1444717851150981'
    }
}

#  Dictionary with corrections of data
RD_BIOBANKS = {
    168144: {
        'name of host institution': 'MRC Centre for Neuromuscular Diseases BioBank London, Dubowitz Neuromuscular Unit'
    },
    168284: {
        'name': 'Muscle Tissue Culture Collection',
        'acronym': 'MTCC'
    },
    168562: {
        'name of host institution': 'Tel-Aviv University'
    },
    44001: {
        'name of host institution': 'Carlos III Health Institute'
    },
    45401: {
        'name': 'Cells, tissues and DNA from patients with neuromuscular diseases',
        'name of host institution': 'Muscle Cell Biology Lab, Neuromuscular Diseases and Neuroimmunology Unit, Ist. Neurologico C. Besta'
    },
    76957: {
        'name': 'Bank for the Diagnosis and Research on Neuromuscular Disorders (NHMGB)'
    },
    77088: {
        'name': 'Genomic and Genetic Disorders Biobank',
        'description': 'The Genomic Disorders Biobank (GDB) located at the “Medical Genetic Unit”, IRCCS Casa Sollievo della Sofferenza, San Giovanni Rotondo, establishes and stores biological samples (DNA, RNA and tissue cell lines) from a number of individuals affected by human rare genomic and genetics diseases. The Biobank includes also molecular and cytogenetics characterization of all samples achieved by cytogenetics and molecular analyses.',
        'url': ['https://www.operapadrepio.it/ggdbbank/index.php'],
        'name of host institution': 'Medical Genetics Unit, Fondazione IRCCS Casa Sollievo della Sofferenza'
    },
    77761: {
        'name of host institution': 'Parkinson Institute, ASST Gaetano Pini CTO (ex Istituti Clinici di Perfezionamento)',
        'url': 'https://www.parkinson.it/biobanca.html'
    }
}

COUNTRIES_CODES = {
    'Belgium': 'BE',
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

NATIONAL_NODES = {
    'Germany': 'DE',
    'Hungary': 'HU',
    'Israel': 'EXT',
    'Italy': 'IT',
    'Malta': 'MT',
    'Slovenia': 'SI',
    'Spain': 'ES',
    'Turkey': 'EXT',
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

RD_NETWORK = 'bbmri-eric:networkID:EXT_RD-Biobanks'

BIOBANKS_SHEET = 'eu_bbmri_eric_biobanks'
RD_BIOBANKS_SHEET = 'eu_bbmri_eric_biobanks_rd'
COLLECTIONS_SHEET = 'eu_bbmri_eric_collections'
RD_COLLECTIONS_SHEET = 'eu_bbmri_eric_collections_rd'
PERSONS_SHEET = 'eu_bbmri_eric_persons'
RD_PERSONS_SHEET = 'eu_bbmri_eric_persons_rd'
ALSO_KNOWN_SHEET = 'eu_bbmri_eric_also_known_in'
RD_ALSO_KNOWN_SHEET = 'eu_bbmri_eric_also_known_in_rd'
DISEASES_SHEET = 'eu_bbmri_eric_disease_types'
COUNTRIES_SHEET = 'eu_bbmri_eric_countries'
NETWORKS_SHEET = 'eu_bbmri_eric_networks'

MISSING_DISEASES_FILE = 'missing_diseases.csv'
MISSING_COUNTRIES_FILE = 'missing_countries.csv'
# MISSING_NETWORKS_FILE = 'missing_networks.csv'


class RDConnectImporter:

    def __init__(self, bb_dir_data, rdc_finder_data, diseases_data, countries_data, networks_data, sc_url, sc_user,
                 sc_pass, add_new_sheets=False):
        self.eric_data = bb_dir_data
        self.add_new_sheets = add_new_sheets
        if add_new_sheets:
            self.eric_data[RD_COLLECTIONS_SHEET] = pd.DataFrame(columns=self.eric_data[COLLECTIONS_SHEET].columns)
            self.eric_data[RD_BIOBANKS_SHEET] = pd.DataFrame(columns=self.eric_data[BIOBANKS_SHEET].columns)
            self.eric_data[RD_PERSONS_SHEET] = pd.DataFrame(columns=self.eric_data[PERSONS_SHEET].columns)
            self.eric_data[RD_ALSO_KNOWN_SHEET] = pd.DataFrame(columns=self.eric_data[ALSO_KNOWN_SHEET].columns)
        self.rdc_finder_data = rdc_finder_data
        self.sc_url = sc_url
        self.sc_user = sc_user
        self.sc_pass = sc_pass
        self.sc_session = Session(self.sc_url)
        self.sc_session.login(self.sc_user, self.sc_pass)

        self.diseases_data = diseases_data
        self.countries_data = countries_data
        self.networks_data = networks_data

        self.missing_diseases = []
        self.missing_countries = []
        # self.check_missing_networks([RD_NETWORK])

    def _concat_cell(self, current_cell, new_values):
        if str(current_cell.values[0]) == 'nan':
            return ','.join(sorted(new_values))
        else:
            new_set = set(current_cell.values[0].split(',')).union(set(new_values))
            return ','.join(sorted(new_set))

    def get_country_code(self, rd_biobank_id, rd_connect_country):
        if rd_biobank_id == 168562:
            return 'IL'
        if rd_biobank_id == 77088:
            return 'IT'
        if rd_biobank_id == 168144:
            return 'UK'
        return COUNTRIES_CODES[rd_connect_country]

    def get_national_nodes_codes(self, rd_biobank_id, rd_connect_country):
        if rd_biobank_id == 168562:
            return ''
        if rd_biobank_id == 77088:
            return 'IT'
        if rd_biobank_id == 168144:
            return 'UK'
        return NATIONAL_NODES[rd_connect_country]

    def biobank_exist(self, rd_biobank_id):
        return rd_biobank_id in BIOBANKS_MAPPING

    def get_biobank_id(self, rd_biobank_id):
        try:
            return BIOBANKS_MAPPING[rd_biobank_id]['biobank']
        except KeyError:
            return f'bbmri-eric:ID:EXT_{rd_biobank_id}'

    def get_collection_id(self, rd_biobank_id):
        try:
            biobank = BIOBANKS_MAPPING[rd_biobank_id]
        except KeyError:
            return f'{self.get_biobank_id(rd_biobank_id)}:collection:MainCollection'
        else:
            try:
                return biobank['collection']
            except KeyError:
                return f'{self.get_biobank_id(rd_biobank_id)}:collection:{rd_biobank_id}'

    def generate_contact_id(self, rd_biobank_id):
        return f'bbmri-eric:contactID:EXT_{rd_biobank_id}'

    def get_material(self, rd_material):
        try:
            return MATERIAL_TYPES[rd_material.upper()]
        except KeyError:
            return None

    def get_donors_and_diseases(self, rd_data):
        num_donors = 0
        diseases = set()
        for d in rd_data['diseases']:
            diseases.update(d['orphacode'])
            diseases.update(icd for icd in d['icd10'] if 'icd10cm' not in icd)
            num_donors += int(d['number'])
        sc_diseases = self.get_diseases_from_sample_catalogue(rd_data['OrganizationID'])
        diseases.update(sc_diseases)
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

    def get_field_value(self, rd_biobank_id, rd_data, field_name):
        try:
            return RD_BIOBANKS[rd_biobank_id][field_name]
        except KeyError:
            return rd_data[field_name]

    def add_biobank_info(self, rd_data):
        rd_biobank_id = rd_data['OrganizationID']
        rd_connect_country = rd_data['address']['country']

        country_code = self.get_country_code(rd_biobank_id, rd_connect_country)
        biobank_id = self.get_biobank_id(rd_biobank_id)

        biobanks_df = self.eric_data[BIOBANKS_SHEET]
        if biobanks_df[biobanks_df.id == biobank_id].empty:  # The biobank is not present in the original data
            new_biobank = pd.DataFrame({
                'id': [biobank_id],
                'pid': [biobank_id],
                'name': [self.get_field_value(rd_biobank_id, rd_data, 'name')],
                'acronym': [self.get_field_value(rd_biobank_id, rd_data['bb_core'], 'acronym')],
                'description': [BeautifulSoup(self.get_field_value(rd_biobank_id, rd_data['bb_core'], 'Description'),
                                              'html.parser').get_text()],
                'url': [self.get_field_value(rd_biobank_id, rd_data, 'url')[0]],
                'location': [None],
                'country': [country_code],
                'head': [None],
                'contact': [self.generate_contact_id(rd_biobank_id)],
                'juridical_person': [
                    self.get_field_value(rd_biobank_id, rd_data['address'], 'name of host institution')],
                'network': [RD_NETWORK],
                'also_known': [None],
                'collections': [self.get_collection_id(rd_biobank_id)],
                'capabilities': [None],
                'quality': [None],
                'collaboration_commercial': [None],
                'collaboration_non_for_profit': [None],
                'national_node': [None],
                'withdrawn': [False]
            }, index=None)
            self.eric_data[BIOBANKS_SHEET] = pd.concat([biobanks_df, new_biobank])
            if self.add_new_sheets:
                self.eric_data[RD_BIOBANKS_SHEET] = pd.concat([self.eric_data[RD_BIOBANKS_SHEET], new_biobank])
        else:
            biobanks_df.loc[biobanks_df.id == biobank_id, 'network'] = \
                self._concat_cell(biobanks_df.loc[biobanks_df.id == biobank_id, 'network'], [RD_NETWORK])
            if self.add_new_sheets:
                self.eric_data[RD_BIOBANKS_SHEET] = pd.concat(
                    [self.eric_data[RD_BIOBANKS_SHEET], biobanks_df.loc[biobanks_df.id == biobank_id]])

    def add_collection_info(self, rd_data):
        rd_biobank_id = rd_data['OrganizationID']
        rd_connect_country = rd_data['address']['country']

        country_code = self.get_country_code(rd_biobank_id, rd_connect_country)
        self.check_missing_country(country_code)
        biobank_id = self.get_biobank_id(rd_biobank_id)
        collection_id = self.get_collection_id(rd_biobank_id)

        materials = rd_data['bb_core']['Biomaterial_Available'] + rd_data['bb_core']['Additional_Biomaterial_available']
        materials = set(self.get_material(m) for m in materials if self.get_material(m) is not None)

        num_of_donors, diseases = self.get_donors_and_diseases(rd_data)

        self.check_missing_diseases(diseases)

        collections_df = self.eric_data[COLLECTIONS_SHEET]
        biobanks_df = self.eric_data[BIOBANKS_SHEET]

        if collections_df[collections_df.id == collection_id].empty:  # If the collection is a new collection
            if rd_biobank_id not in BIOBANKS_MAPPING.keys():
                # If also the biobank was not present it means it is a completely new collection/biobank
                # so the name will be given to the Biobank and the collection will be given the generic name
                name = 'Main Collection'
                acronym = ''
                description = ''
                contact_id = self.generate_contact_id(rd_biobank_id)
                biobank_label = self.get_field_value(rd_biobank_id, rd_data, 'name')
                combined_network = RD_NETWORK
                #  In this case it also creates the biobank
                self.add_biobank_info(rd_data)
                self.add_contact_info(rd_data)
            else:
                # If the biobank is present, it means this is a new collection of an old biobank
                # so the collection will be given the name of the finder biobank
                name = self.get_field_value(rd_biobank_id, rd_data, 'name')
                acronym = self.get_field_value(rd_biobank_id, rd_data['bb_core'], 'acronym')
                description = BeautifulSoup(self.get_field_value(rd_biobank_id, rd_data['bb_core'], 'Description'),
                                            'html.parser').get_text()
                contact_id = biobanks_df.loc[biobanks_df.id == biobank_id, 'contact'].values[0]
                biobank_label = biobanks_df.loc[biobanks_df.id == biobank_id, 'name'].values[0]
                combined_network = self._concat_cell(
                    biobanks_df.loc[biobanks_df.id == biobank_id, 'network'], [RD_NETWORK])

                # it appends the collection to the list of collections of the biobank
                biobanks_df.loc[biobanks_df.id == biobank_id, 'collections'] = \
                    self._concat_cell(biobanks_df.loc[biobanks_df.id == biobank_id, 'collections'], [collection_id])
                if self.add_new_sheets:
                    self.eric_data[RD_BIOBANKS_SHEET] = pd.concat(
                        [self.eric_data[RD_BIOBANKS_SHEET], biobanks_df.loc[biobanks_df.id == biobank_id]])
                self.add_contact_info(rd_data)

            new_collection = pd.DataFrame({
                'id': [collection_id],
                'name': [name], 'acronym': [acronym], 'description': [description], 'url': [None], 'location': [None],
                'country': [country_code], 'head': [None],
                'contact': [contact_id],
                'withdrawn': [False], 'national_node': [None],
                'parent_collection': [None], 'sub_collections': [None],
                'biobank': [biobank_id], 'biobank_label': [biobank_label], 'network': [RD_NETWORK],
                'combined_network': [combined_network], 'also_known': [f'rdconnect:{rd_biobank_id}'],
                'type': ['RD'],  # statically add only rare disease type
                'data_categories': ['BIOLOGICAL_SAMPLES,OTHER'],
                'order_of_magnitude': ['0'],  # TODO: check if it can be retrieved in the samples catalogu]e
                'size': [None],  # TODO: check if it can be retrieved in the samples catalogue
                'categories': [None], 'timestamp': [None], 'quality': [None], 'combined_quality': [None],
                'number_of_donors': [num_of_donors],
                'order_of_magnitude_donors': [int(math.log10(max(1, num_of_donors)))],
                'sex': [None],  # TODO: check if it can be retrieved in the samples catalogue
                'diagnosis_available': [','.join(sorted(diseases))], 'age_low': [None],
                'age_high': [None], 'age_unit': [None], 'materials': [','.join(sorted(materials))],
                'storage_temperatures': [None], 'body_part_examined': [None], 'imaging_modality': [None],
                'image_dataset_type': [None], 'collaboration_commercial': [None],
                'collaboration_non_for_profit': [None],
                'data_use': [None], 'commercial_use': [None], 'access_fee': [None], 'access_joint_project': [None],
                'access_description': [None], 'access_uri': [None], 'sop': [None]
            })
            self.eric_data[COLLECTIONS_SHEET] = pd.concat([collections_df, new_collection])
            if self.add_new_sheets:
                self.eric_data[RD_COLLECTIONS_SHEET] = pd.concat([self.eric_data[RD_COLLECTIONS_SHEET], new_collection])
        else:
            # If the collection was already present in the Directory, it updates some fields
            collections_df.loc[collections_df.id == collection_id, ['also_known', 'network', 'combined_network',
                                                                    'diagnosis_available']] = \
                [self._concat_cell(collections_df.loc[collections_df.id == collection_id, 'also_known'],
                                   [f'rdconnect:{rd_biobank_id}']),
                 self._concat_cell(collections_df.loc[collections_df.id == collection_id, 'network'], [RD_NETWORK]),
                 self._concat_cell(collections_df.loc[collections_df.id == collection_id, 'combined_network'],
                                   [RD_NETWORK]),
                 self._concat_cell(collections_df.loc[collections_df.id == collection_id, 'diagnosis_available'],
                                   sorted(diseases))]

            if self.add_new_sheets:
                # Add the rows
                self.eric_data[RD_COLLECTIONS_SHEET] = pd.concat(
                    [self.eric_data[RD_COLLECTIONS_SHEET], collections_df.loc[collections_df.id == collection_id]])
                self.eric_data[RD_BIOBANKS_SHEET] = pd.concat(
                    [self.eric_data[RD_BIOBANKS_SHEET], biobanks_df.loc[biobanks_df.id == biobank_id]])

    def add_contact_info(self, rd_data):
        rd_biobank_id = rd_data['OrganizationID']
        rd_connect_country = rd_data['address']['country']
        country_code = self.get_country_code(rd_biobank_id, rd_connect_country)
        contact_id = self.generate_contact_id(rd_biobank_id),

        df = self.eric_data[PERSONS_SHEET]
        if df[df.id == contact_id].empty:
            biobank_id = self.get_biobank_id(rd_biobank_id)
            collection_id = self.get_collection_id(rd_biobank_id)
            new_contact = pd.DataFrame({
                'id': [self.generate_contact_id(rd_biobank_id)],
                'title_before_name': [None],
                'first_name': [rd_data['main contact']['first name']],
                'last_name': [rd_data['main contact']['last name']],
                'title_after_name': [None],
                'email': [rd_data['main contact']['email']],
                'phone': [None],
                'address': [None],
                'zip': [None],
                'city': [None],
                'country': [country_code],
                'role': [None],
                'biobanks': [biobank_id],
                'collections': [collection_id],
                'networks': [None],
                'national_node': [None],
                'withdrawn': [None]
            })
            self.eric_data[PERSONS_SHEET] = pd.concat([df, new_contact])
            if self.add_new_sheets:
                self.eric_data[RD_PERSONS_SHEET] = pd.concat([self.eric_data[RD_PERSONS_SHEET], new_contact])

    def add_also_known_in_info(self, rd_data):
        rd_biobank_id = rd_data['OrganizationID']

        df = self.eric_data[ALSO_KNOWN_SHEET]
        also_known_id = f'rdconnect:{rd_biobank_id}'
        new_also_known = pd.DataFrame({
            'id': [also_known_id],
            'name_system': ['RD Connect'],
            'pid': [rd_biobank_id],
            'url': [None],
            'national_node': [self.get_national_nodes_codes(rd_biobank_id, rd_data['address']['country'])],
            'withdrawn': [None]
        })
        self.eric_data[ALSO_KNOWN_SHEET] = pd.concat([df, new_also_known])
        if self.add_new_sheets:
            self.eric_data[RD_ALSO_KNOWN_SHEET] = pd.concat([self.eric_data[RD_ALSO_KNOWN_SHEET], new_also_known])

    def check_missing_diseases(self, diseases):
        disease_df = self.diseases_data
        for d in diseases:
            if self.eric_data[DISEASES_SHEET][self.eric_data[DISEASES_SHEET].id == d].empty:
                if disease_df is not None and not disease_df[disease_df.id == d].empty:
                    self.eric_data[DISEASES_SHEET] = pd.concat(
                        [self.eric_data[DISEASES_SHEET], disease_df[disease_df.id == d]])
                else:
                    self.missing_diseases.append(d)

    def check_missing_country(self, country):
        country_df = self.countries_data
        if self.eric_data[COUNTRIES_SHEET][self.eric_data[COUNTRIES_SHEET].id == country].empty:
            if country_df is not None and not country_df[country_df.id == country].empty:
                self.eric_data[COUNTRIES_SHEET] = pd.concat(
                    [self.eric_data[COUNTRIES_SHEET], country_df[country_df.id == country]])
            else:
                self.missing_countries.append(country)

    # def check_missing_networks(self, networks):
    #     networks_df = self.networks_data
    #     for n in networks:
    #         if self.eric_data[NETWORKS_SHEET][self.eric_data[NETWORKS_SHEET].id == n].empty:
    #             if networks_df is not None and not networks_df[networks_df.id == n].empty:
    #                 self.eric_data[NETWORKS_SHEET] = pd.concat(
    #                     [self.eric_data[NETWORKS_SHEET], networks_df[networks_df.id == n]])

    def run(self):
        for b in self.rdc_finder_data:
            print("Converting biobank: ", b['OrganizationID'])
            print("Getting collection data: ", b['OrganizationID'])
            self.add_collection_info(b)
            # print("Getting person data: ", b['OrganizationID'])
            # self.add_contact_info(b)
            print("Getting also known data: ", b['OrganizationID'])
            self.add_also_known_in_info(b)
            print()
        if len(self.missing_diseases) > 0:
            print("Some disease codes could not be added. Add the following codes manually")
            print(sorted(set(self.missing_diseases)))
        return self.eric_data


def write_excel(data, output_file):
    print("Saving excel")
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
                        help='Output EMX file containing updated data', default='bbmri-directory-with-rd-biobanks.xlsx')
    parser.add_argument('--sc-url', '-H', dest='sc_url', type=str, required=True,
                        help='Base URL of the RD Connect Samples Catalogue')
    parser.add_argument('--sc-user', '-u', dest='sc_user', type=str, required=True,
                        help='Samples Catalogue user name')
    parser.add_argument('--sc-password', '-p', dest='sc_pwd', type=str, required=True,
                        help='Samples Catalogue user password')


    args = parser.parse_args()

    print('loading source')
    with open(args.rd_connect_file) as f:
        rd_connect_data = json.load(f)

    with open(MISSING_DISEASES_FILE) as f:
        diseases = pd.read_csv(f, sep=';', header=0, index_col=False)

    with open(MISSING_COUNTRIES_FILE) as f:
        countries = pd.read_csv(f, sep=';', header=0, index_col=False)

    networks = []
    rdc_biobanks = sorted(rd_connect_data['allData'], key=lambda b: b['OrganizationID'])
    rdc_biobanks = filter(lambda b: b['OrganizationID'] in ITALIAN_BIOBANKS, rdc_biobanks)
    d_emx = pd.read_excel(args.directory_file, sheet_name=None, engine='openpyxl')
    importer = RDConnectImporter(d_emx, rdc_biobanks, diseases, countries, networks, args.sc_url, args.sc_user,
                                 args.sc_pwd, add_new_sheets=True)
    res = importer.run()
    write_excel(res, args.output_file)
