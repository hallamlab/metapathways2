"""This module provides classes for dealing with PathoLogic files."""

import re

from __init__ import SequenceRecordSerializer

PRODUCT_TYPE_MAPPINGS = {
    'pseudo': 'PSEUDO',
    'tRNA': 'TRNA',
    'rRNA': 'RRNA',
    'misc_RNA': 'MISC-RNA',
}

DATABASE_MAPPINGS = {
    'InterPro': 'UNIPROT',
    'GI': 'GI',
    'GeneID': 'GeneID',
    'GO': 'GO'
}

untagged_i = 0

def split_ec_numbers(ec_number):
    ec_numbers = ec_number.split(" ")
    results = []
    for ec_number in ec_numbers:
       if len(ec_number.strip()) > 0:
          results.append(ec_number)
    return results

def feature_to_string(feature):
    """Converts a specified feature to a PathoLogic entry.
    
    Several entries may be created if the sequence of the entry is composed of
    joins.
    """

    global untagged_i

    # Initialization
    details = {
        'product_type': None,
        'references': [],
        'go_terms': [],
        'ec_numbers': []
    }

    # Determine the product type

    details['product_type'] = (PRODUCT_TYPE_MAPPINGS[feature.type] if
                               feature.type in PRODUCT_TYPE_MAPPINGS else 'P')

    # ID and name

    if feature.locus_tag:
        details['id'] = feature.locus_tag
        details['name'] = feature.locus_tag
    else:
        details['id'] = 'Untagged_%04d' % untagged_i
        details['name'] = 'Untagged_%04d' % untagged_i
        untagged_i += 1

    # EC numbers

    details['ec_numbers'] = set(feature.ec_numbers)

    # Qualified name

    if feature.product:
        details['function'] = feature.product
    elif feature.functions:
        details['function'] = ''.join(feature.functions)
    else:
        details['function'] = 'ORF'

    details['product'] = feature.product

    if feature.notes:
        for field in '\n'.join(feature.notes).split('\n'):

            matches = re.match(r'^KEGG:\s+[a-z]{3}:\S+\s+(.*)', field)
            if matches:
                details['function'] = matches.group(1)

            matches = re.match(r'^PFAM:\s+(.*)', field)
            if matches:
                details['function'] = matches.group(1)

    # Coordinates

    details['coordinates'] = [(str(c.begin), str(c.end)) for c in
                              feature.location.coordinates]

    # Database cross-references

    details['references'] = []

    if feature.db_xrefs:
        for xref in feature.db_xrefs:
            database, reference = xref.split(':')
            if database in DATABASE_MAPPINGS:
                details['references'].append('%s:%s' %
                                             (DATABASE_MAPPINGS[database],
                                              reference))

    # GO terms

    details['go_terms'] = ([go.split(':')[1] for go in feature.db_xrefs
                            if go.startswith('GO:')])

    # Clean up the function to remove parentheses and newlines

    if 'function' in details:
        details['function'] = re.sub(r'\(.*\)', '', details['function'])
        details['function'] = re.sub(r'\n', ' ', details['function'])

    # Create the entries...

    def create_entry(coord_tuple=None, i=None):

        entry_lines = []

        entry_id = '%s_%i' % (details['id'], i) if i else details['id']

        if 'id' in details:
            entry_lines.append('ID\t%s' % entry_id)

        if 'name' in details:
            entry_lines.append('NAME\t%s' % details['name'])

        if 'function' in details:
            entry_lines.append('FUNCTION\t%s' % details['function'])

        if  re.search("hypothetical protein",details['function']):
            return("")

        if coord_tuple:
            entry_lines.append('STARTBASE\t%s' % coord_tuple[0])
            entry_lines.append('ENDBASE\t%s' % coord_tuple[1])

        if 'product_type' in details:
            entry_lines.append('PRODUCT-TYPE\t%s' % details['product_type'])

        if 'references' in details:
            for reference in details['references']:
                entry_lines.append('DBLINK\t%s' % reference)

        if 'go_terms' in details:
            for go_term in details['go_terms']:
                entry_lines.append('GO\t%s' % go_term)

        if 'ec_numbers' in details:
            for ec_number in details['ec_numbers']:
                ec_numbers = split_ec_numbers(ec_number)
                for ec_number in ec_numbers:
                   entry_lines.append('EC\t%s' % ec_number)

        return '\n'.join(entry_lines)

    entries = []

    if 'coordinates' in details:
        if len(details['coordinates']) > 1:
            for i, coordinate in enumerate(details['coordinates']):
                entry =create_entry(coordinate, i + 1)
                if len(entry) > 0:
                   entries.append(create_entry(coordinate, i + 1))
        else:
            entry = create_entry(details['coordinates'][0])
            if len(entry) > 0:
                 entries.append(entry)

    else:
        entry = create_entry()
        if len(entry) > 0:
             entries.append(entry)


    string_of_entries = ""
    for entry  in entries:
       if len(string_of_entries) > 0 :
           string_of_entries =  string_of_entries + '\n//\n' + entry
       else:
           string_of_entries = entry


    return  string_of_entries


class PathoLogicRecordSerializer(SequenceRecordSerializer):
    """Serializes a record to PathoLogic format."""
    def serialize(self, record):
        filtered_features = (feature for feature in record.features if
                             feature.type in ('CDS', 'tRNA', 'rRNA', 'mRNA',
                                              'misc_RNA'))
        results = [] # feature_to_string(f) for f in filtered_features]

        for f in filtered_features:
           #print feature_to_string(f)
           isEmpty = True 
           for gene in feature_to_string(f):
              if  len(gene) > 0:
                isEmpty = False

           if not isEmpty:
              results.append(feature_to_string(f))

        return '\n//\n'.join(results)





