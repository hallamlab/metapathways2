"""This module defines classes for working with FASTA files."""

import textwrap

from __init__ import (SequenceRecordParser, SequenceRecordSerializer,
                      SequenceRecord, SequenceFactory, Feature)
import location

class FastaRecordParser(SequenceRecordParser):
    """Parses a FASTA record from a string."""

class FastaRecordSerializer(SequenceRecordSerializer):
    """Serializes a record to a string in FASTA format.
    
    You may customize the format of the description by specifying the
    *description_format*. Possible substitutions include:

    =================   ================
    Parameter           What it is
    =================   ================
    ``%(gi_number)s``   GI number
    ``%(accession)s``   Accession number
    ``%(locus)s``       Locus
    =================   ================

    .. note::

        If you are specifying a custom description format, avoid formats that
        result in excessively long descriptions.
    """
    SEQUENCE_WRAPPER = textwrap.TextWrapper(width=79)

    def __init__(self,
                 description_format='gi|%(gi_number)s|gb|%(accession)s|%(locus)s'):
        self.description_format = description_format

    def serialize(self, record):
        description = self.description_format % {'gi_number': None,
                                                 'accession': record.accession,
                                                 'locus': record.locus}
        sequence = '\n'.join(self.SEQUENCE_WRAPPER.wrap(str(record.sequence)))
        return '>%s\n%s' % (description, sequence)
