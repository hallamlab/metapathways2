"""This module defines classes for working with sequences and records."""

import datetime
import re
import location

# Sequences

class Sequence(object):
    """Defines a biological sequence."""
    def __init__(self, sequence):
        self.sequence = sequence

    def __iter__(self):
        return iter(self.sequence)

    def __str__(self):
        return self.sequence

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        truncated = (self.sequence.replace(self.sequence[10:-10], '...')
                     if len(self.sequence) > 20 else self.sequence)
        return '<%s %s length=%d>' % (self.__class__.__name__, truncated,
                                      len(self.sequence))

    def __getslice__(self, start, end):
        return self.sequence[start:end]

    def __add__(self, other):
        if type(self) != type(other):
            raise TypeError()
        return type(self)(self.sequence + other.sequence)



class SequenceFactory(object):
    def get_sequence(self, sequence, sequence_type):
        sequence_type = sequence_type.upper()
        if sequence_type == 'DNA':
            return DNASequence(sequence)
        elif sequence_type == 'RNA':
            return RNASequence(sequence)
        elif sequence_type == 'AA':
            return AminoAcidSequence(sequence)
        else:
            raise ValueError('Unknown sequence type %s' % repr(sequence_type))



class NucleotideSequence(Sequence):
    """Defines a nucleotide sequence."""
    pairings = {}
    type = None

    def __init__(self, sequence):
        self.sequence = sequence.upper()

    def __get_reverse_strand(self):
        return self.__class__(''.join(self.pairings[base] if base in
                                      self.pairings else '?' for base in
                                      reversed(self.sequence)))
    #: The reverse strand, in the 5' -> 3' direction.
    reverse_strand = property(__get_reverse_strand)



class RNASequence(NucleotideSequence):
    """Defines an RNA sequence."""
    type = 'RNA'
    pairings = {
        'A': 'U',
        'U': 'A',
        'G': 'C',
        'C': 'G',
        'N': 'N'
    }



class DNASequence(NucleotideSequence):
    """Defines a DNA sequence."""
    type = 'DNA'
    pairings = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'N': 'N'
    }



class AminoAcidSequence(Sequence):
    """Defines an amino acid sequence."""



# Sequence records

class SequenceRecord(object):
    """Defines a record for a sequence."""
    # TODO: Add GI number
    def __init__(self, locus=None, accession=None, division=None,
                 molecule_type=None, date=None, definition=None,  version=None,
                 keywords=None, source=None, source_organisms=None,
                 references=None, comments=None, features=None, sequence=None):
        self.locus = locus
        self.accession = accession if accession else []
        self.division = division
        self.molecule_type = molecule_type
        self.date = date if date else datetime.datetime.now()
        self.definition = definition
        self.version = version
        self.keywords = keywords if keywords else []
        self.source = source
        self.source_organisms = source_organisms if source_organisms else []
        self.references = references if references else []
        self.comments = comments if comments else []

        self.features = features if features else []
        self.sequence = sequence


    def __set_references(self, value):
        self.references = value

    references_ = property(fset= __set_references)

    def __get_translations(self):
        for feature in self.features:
            qualifiers = feature.qualifiers
            if ('locus_tag' in qualifiers and 'translation' in qualifiers
                and qualifiers['locus_tag'] and qualifiers['translation']):
                
                # Generate a sequence and source feature
                print qualifiers['translation']
                sequence = AminoAcidSequence(''.join(qualifiers['translation']))
                feature = Feature('source',
                                  location.Location([location.Coordinates(location.Coordinate(1),
                                                                          location.Coordinate(len(sequence)))]))

                # Create the record
                record = SequenceRecord(locus=''.join(qualifiers['locus_tag']),
                                        accession=self.accession,
                                        division=self.division,
                                        molecule_type=self.molecule_type,
                                        date=self.date,
                                        definition=(''.join(qualifiers['product'])
                                                    if 'product' in qualifiers
                                                    and qualifiers['product']
                                                    else None),
                                        version=self.version,
                                        keywords=[],
                                        source=self.source,
                                        source_organisms=self.source_organisms,
                                        references=self.references,
                                        comments=self.comments,
                                        features=[feature],
                                        sequence=sequence)
                yield record
        return
    #: A list of :class:`SequenceRecord` objects containing translation
    #: sequences.
    translations = property(__get_translations)

    def __get_protein_encoding_sequences(self):
        # TODO: If we ever need to, we can adjust this to figure out the first
        # triplet based on the translations
        for feature in self.features:
            qualifiers = feature.qualifiers
            if ('locus_tag' in qualifiers and 'translation' in qualifiers
                and qualifiers['locus_tag'] and qualifiers['translation']):
                
                # Generate a DNA sequence

                # ...by collecting a list of coordinates
                # TODO: Deal with joins properly
                sequence_parts = []
                for coordinate in feature.location.coordinates:
                    begin = int(coordinate.begin)
                    end = int(coordinate.end)
                    if begin < end:
                        sequence_parts.append(self.sequence[begin - 1:end])
                    else:
                        sequence_parts.append(self.sequence[end - 1:begin])

                sequence = (DNASequence(''.join(sequence_parts)))
                feature = Feature('source',
                                  location.Location([location.Coordinates(location.Coordinate(1),
                                                                          location.Coordinate(len(sequence)))]))
                
                # Create the record
                record = SequenceRecord(locus=''.join(qualifiers['locus_tag']),
                                        accession=self.accession,
                                        division=self.division,
                                        molecule_type=self.molecule_type,
                                        date=self.date,
                                        definition=(''.join(qualifiers['product'])
                                                    if 'product' in qualifiers
                                                    and qualifiers['product']
                                                    else None),
                                        version=self.version,
                                        keywords=[],
                                        source=self.source,
                                        source_organisms=self.source_organisms,
                                        references=self.references,
                                        comments=self.comments,
                                        features=[feature],
                                        sequence=sequence)
                yield record
    protein_encoding_sequences = property(__get_protein_encoding_sequences)



    def __repr__(self):
        return (('<SequenceRecord locus=%s accession=%s division=%s '
                 'molecule_type=%s date=%s definition=%s version=%s '
                 'keywords=%s source=%s source_organisms=%s references=%s '
                 'features=%s sequence=%s>') % (self.locus, self.accession,
                 self.molecule_type, self.division, repr(self.date),
                 repr(self.definition), repr(self.version),
                 repr(self.keywords), repr(self.source),
                 repr(self.source_organisms), repr(self.references),
                 repr(self.features), repr(self.sequence)))



class Attribute(object):
    def __init__(self, name, value, subattributes=None):
        self.name = name
        self.value = value
        self.subattributes = subattributes if subattributes else []

    def __repr__(self):
        return '<Attribute %s=%s subattributes=%s>' % (self.name,
                                                       repr(self.value),
                                                       repr(self.subattributes))



class Feature(object):
    COMPLEMENT_PATTERN = re.compile(r'complement', re.I)
    COORDINATES_PATTERN = re.compile(r'(\d+)[.]*(\d+)', re.I)

    def __init__(self, type, location, qualifiers=None):
        self.type = type
        self.location = location
        self.qualifiers = qualifiers if qualifiers else {}

    def __get_locus_tag(self):
        if 'locus_tag' in self.qualifiers and self.qualifiers['locus_tag']:
            return self.qualifiers['locus_tag'][0]
        return None

    def __set_locus_tag(self, value):
        self.qualifiers['locus_tag'] = [value]

    locus_tag = property(__get_locus_tag, __set_locus_tag)

    def __get_product(self):
        if 'product' in self.qualifiers and self.qualifiers['product']:
            return self.qualifiers['product'][0]
        return None

    def __set_product(self, product):
        if not isinstance(product, (str, unicode)) and product is not None:
            raise TypeError('The product %s is not a string' % repr(product))
        self.qualifiers['product'] = [product] if product else []

    product = property(__get_product, __set_product)

    def __get_db_xrefs(self):
        if 'db_xref' not in self.qualifiers:
            self.qualifiers['db_xref'] = set()
            #self.qualifiers['db_xref'] = []
        self.qualifiers['db_xref'] = set(self.qualifiers['db_xref'])
        return self.qualifiers['db_xref']

    def __set_db_xrefs(self, xrefs):
        self.qualifiers['db_xref'] = set(xrefs)

    db_xrefs = property(__get_db_xrefs, __set_db_xrefs)

    def __get_ec_numbers(self):
        if 'EC_number' not in self.qualifiers:
            self.qualifiers['EC_number'] = set()
        self.qualifiers['EC_number'] = set(self.qualifiers['EC_number'])
        return self.qualifiers['EC_number']

    def __set_ec_numbers(self, numbers):
        self.qualifiers['EC_number'] = list(numbers)

    ec_numbers = property(__get_ec_numbers, __set_ec_numbers)

    def __get_notes(self):
        if 'note' not in self.qualifiers:
            self.qualifiers['note'] = []
        return self.qualifiers['note']

    notes = property(__get_notes)

    def __get_translation(self):
        if 'translation' in self.qualifiers and self.qualifiers['translation']:
            return self.qualifiers['translation'][0]
        return None

    def __set_translation(self, txn):
        if not isinstance(txn, Sequence):
            raise TypeError('The translation %s is not a Sequence' % repr(txn))
        self.qualifiers['translation'] = [txn]

    translation = property(__get_translation, __set_translation)

    def __get_organism(self):
        if 'organism' in self.qualifiers and self.qualifiers['organism']:
            return self.qualifiers['organism'][0]
        return None
    
    def __set_organism(self, value):
        if not isinstance(value, (str, unicode)):
            raise TypeError('The organism %s is not a string' % repr(value))
        self.qualifiers['organism'] = [value]

    organism = property(__get_organism, __set_organism)

    def __get_strain(self):
        if 'strain' in self.qualifiers and self.qualifiers['strain']:
            return self.qualifiers['strain'][0]
        return None
    
    def __set_strain(self, value):
        if not isinstance(value, (str, unicode)):
            raise TypeError('The strain %s is not a string' % repr(value))
        self.qualifiers['strain'] = [value]

    strain = property(__get_strain, __set_strain)

    def __get_genes(self):
        if 'gene' not in self.qualifiers:
            self.qualifiers['gene'] = set()
        else:
            self.qualifiers['gene'] = set(self.qualifiers['gene'])
        return self.qualifiers['gene']

    genes = property(__get_genes)

    def __get_functions(self):
        if 'function' not in self.qualifiers:
            self.qualifiers['function'] = []
        return self.qualifiers['function']

    def __set_functions(self, value):
        self.qualifiers['function'] = value
    
    functions = property(__get_functions)

    def __get_coordinates(self):
        strand = '+'
        #m = re.search(r'complement',self.location)
        m = self.COORDINATES_PATTERN.search(str(self.location))
        if  self.COMPLEMENT_PATTERN.search(str(self.location)):
           return (m.group(1), m.group(2), '-')
        else:
           return (m.group(1), m.group(2), '+')

    coordinates = property(__get_coordinates)

    def __repr__(self):
        return '<Feature %s at %s qualifiers=%s>' % (self.type, self.location,
                                                     repr(self.qualifiers))



class SequenceRecordParser(object):
    """Parses a sequence record from a string."""
    def __init__(self, contents):
        self.contents = contents



class SequenceRecordSerializer(object):
    """Writes a sequence record to a string."""
