"""This module defines classes for working with InterProScan output."""

import datetime
import re

class InterProScanParser(object):
    """Parses InterProScan output."""
    GO_TERM_REGEXP = re.compile(r'^(.*?): (.*?) \((GO:\d{7})\)$')

    def __init__(self, contents):
        self.contents = contents

    def __get_rows(self):
        input_rows = self.contents.split('\n')
        for input_row in input_rows:
            if input_row:
                (locus_tag, uid, unknown_1, database, accession_number,
                 classification, start_coordinate, end_coordinate, e_value,
                 true_positive, date, interpro_id, protein_name,
                 go_terms) = input_row.split('\t')

                # Massage the input
                unknown_1 = int(unknown_1)
                coordinates = (int(start_coordinate), int(end_coordinate))
                e_value = float(e_value)
                true_positive = true_positive == 'T'
                date = datetime.datetime.strptime(date.upper(), '%d-%b-%Y')
                interpro_id = (interpro_id if interpro_id
                               and interpro_id != 'NULL' else None)
                protein_name = (protein_name if protein_name
                                and protein_name != 'NULL' else None)

                split_go_terms = set()
                for f in go_terms.split(', '):
                    matches = self.GO_TERM_REGEXP.match(f)
                    if matches:
                        split_go_terms.add(matches.groups())

                yield {
                    'locus_tag': locus_tag,
                    'uid': uid,
                    'unknown_1': unknown_1,
                    'database': database,
                    'accession_number': accession_number,
                    'classification': classification,
                    'coordinates': coordinates,
                    'expect_value': e_value,
                    'true_positive': true_positive,
                    'date': date,
                    'interpro_id': interpro_id,
                    'protein_name': protein_name,
                    'go_terms': split_go_terms
                }
    rows = property(__get_rows)
