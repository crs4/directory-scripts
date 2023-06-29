from molgenis import client
from molgenis.client import MolgenisRequestError

ALSO_KNOWN_ENTITY = 'eu_bbmri_eric_also_known_in'
BIOBANK_ENTITY = 'eu_bbmri_eric_biobanks'
COLLECTION_ENTITY = 'eu_bbmri_eric_collections'
CONTACT_ENTITY = 'eu_bbmri_eric_persons'
NETWORK_ENTITY = 'eu_bbmri_eric_networks'

QUERY_ARGS = {
    ALSO_KNOWN_ENTITY: {},
    BIOBANK_ENTITY: {},
    COLLECTION_ENTITY: {'sort_column': 'parent_collection', 'sort_order': 'desc'},
    CONTACT_ENTITY: {},
    NETWORK_ENTITY: {}
}

def add_records(session, entity, records):
    created_records = []
    for i in range(0, len(records), 1000):
        try:
            created_records.extend(session.add_all(entity, records[i:i + 1000]))
        except MolgenisRequestError as ex:
            print("Error adding records")
            print(ex)
    print("Added {} record(s) of type {}".format(len(created_records), entity))


def delete_records(session, entity, records_ids):
    removed_records = []
    for i in range(0, len(records_ids), 1000):
        try:
            removed_records.extend(
                session.delete_list(entity, records_ids[i:i + 1000]))
        except MolgenisRequestError as ex:
            print("Error removing records")
            print(ex)
    print(f"Removed {len(removed_records)} of type {entity}")

def get_attributes_metadata(source_session, entity_type):
    entity_metadata = source_session.get_entity_meta_data(entity_type)
    attributes = []
    # it generates correct description of attributes
    for a in entity_metadata['attributes'].values():
        if a['fieldType'] == 'COMPOUND':  # for compound attributes we need to get metadata of the children attributes
            for ca in a['attributes']:
                new_attr = source_session.get_attribute_meta_data(entity_type, ca['href'].split('/')[-1])
                # the refEntity key of reference types' attributes doesn't have the idAttribute needed, so we need to get it
                if new_attr['fieldType'] in ('XREF', 'CATEGORICAL', 'MREF', 'CATEGORICAL_MREF'):
                    ref_entity = source_session.get_entity_meta_data(
                        new_attr['refEntity']['hrefCollection'].split('/')[-1])
                    new_attr['refEntity']['idAttribute'] = ref_entity['idAttribute']
                attributes.append(new_attr)
        else:
            attributes.append(a)
    return attributes

def sync_record(source_session, dest_session, entity_type):
    attributes = get_attributes_metadata(source_session, entity_type)
    source_records = source_session.get(entity_type, **QUERY_ARGS[entity_type])
    dest_records = []
    for sr in source_records:
        dest_record = {}
        for a in attributes:
            attr_name = a['name']
            attr_type = a['fieldType']
            if attr_name not in sr or attr_type == 'ONE_TO_MANY':
                # skip ONE_TO_MANY since it means it is a readonly attr which is automatically computed
                pass
            elif attr_type in (
                    'STRING', 'TEXT', 'BOOL', 'HYPERLINK', 'DECIMAL', 'LONG', 'INT', 'EMAIL', 'DATE', 'DATE_TIME'):
                dest_record[attr_name] = sr[attr_name]
            elif attr_type in ('XREF', 'CATEGORICAL'):
                dest_record[attr_name] = sr[attr_name][a['refEntity']['idAttribute']]
            elif attr_type in ('MREF', 'CATEGORICAL_MREF'):
                dest_record[attr_name] = [v[a['refEntity']['idAttribute']] for v in sr[attr_name]]
            else:
                print(attr_name)
        dest_records.append(dest_record)

    add_records(dest_session, entity_type, dest_records)

def sync(source_session, dest_session, **kwargs):
    entities = (CONTACT_ENTITY, ALSO_KNOWN_ENTITY, NETWORK_ENTITY, BIOBANK_ENTITY, COLLECTION_ENTITY)

    for entity_type in entities[::-1]:
        old_records = dest_session.get(entity_type)
        delete_records(dest_session, entity_type, [r['id'] for r in old_records])

    for entity_type in entities:
        sync_record(source_session, dest_session, entity_type)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--source-url', '-s')
    parser.add_argument('--dest-url', '-d')
    args = parser.parse_args()
    source = client.Session(args.source_url)
    dest = client.Session(args.dest_url) # if needed, add dest.login or molgenisToken if inside molgenis

    sync(source, dest)
