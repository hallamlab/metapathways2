"""This module defines classes for working with GenBank records."""

import re
import datetime
import textwrap
import sys

from __init__ import (SequenceRecordParser, SequenceRecordSerializer,
                      SequenceRecord, SequenceFactory, Feature)
import location

class GenBankRecordParser(SequenceRecordParser):
    """Parses a GenBank record from a string or file."""
    ATTRIBUTE_REGEXP = re.compile(r'^[A-Z]+')
    SUBATTRIBUTE_REGEXP = re.compile(r'^\s{2,}[A-Z]+')
    FEATURE_REGEXP = re.compile(r'^\s{5}[A-Za-z]+')
    QUALIFIER_REGEXP = re.compile(r'^\s{21}/[A-Za-z_]+')

    FUZZY_LOCATION_REGEXP = re.compile(r'^(\D+)\d+$')

    LOCUS_REGEXP = re.compile(r'^(\S+)\s+\d+ bp\s+(\S+)\s+(.+?)\s+(\d{2}-[A-Z]{3}-\d{4})$')
    count = 1 

    def __init__(self, contents):
        try:
            self.contents = contents.read()
        except AttributeError:
            self.contents = contents
        self.entries = []

    def __iter__(self):

        self.entries = [e for e in self.contents.split('\n//') if e.strip()]
        return self

    def next(self):
        # TODO: Clean this up nicely :)
      #  if self.count % 10000 == 0 :
      #     print "Next " + str(self.count)
      #  self.count = self.count + 1;
        if self.entries:
            lines = self.entries.pop(0).strip().split('\n')
            # Ok, we now have a single genbank entry starting with LOCUS and ending with //
            # one line per element in the list called 'lines'
            
            sequence_type = None
            sequence_record = SequenceRecord()

            # For adding /locus_tag=$locus_$ltnum for CDSes without them:
            ltnum = 1

            # Parse out attributes (the section before FEATURES)
            while (lines and not lines[0].startswith('FEATURES')):
                line = lines.pop(0)  # Get the first line, remove it from the big list

                attribute = line[0:11].strip()  # The first 12 characters, spaces removed
                value = [line[12:]]   # The rest of the line

                # Append lines to 'value' while the lines start with 11 spaces
                while lines and lines[0].startswith(' ' * 11):
                    value.append(lines.pop(0)[12:])

                # Join lines with newlines between them:
                value = '\n'.join(value)
                
                # Do something about the attribute, depending on what it is
                if attribute == 'LOCUS':

                    fields = re.split(r'\s+', value.strip())
                    
                    #for index in range(len(fields)):
                    #    print >>sys.stderr, 'fields[%d]: %s' % ( index,fields[index] )

                    # Two examples seen:
                    #              0                    1       2     3       4        5   6
                    # "LOCUS       U00096               4639675 bp    DNA     circular BCT 14-MAY-2010"
                    #
                    #              0                    1     2       3         4
                    # "LOCUS       amd_scaffold_35      46630 bp      DNA       19-APR-2011 "
                    #  (note trailing space)

                    if fields[2] != 'bp':
                        raise Exception("'bp' expected, got '%s' instead.  Line:\n%s" % (fields[2],value))

                    sequence_record.locus = fields[0]
                    # sequence_record.size = fields[1]  <-- calculated later
                    # 'bp' = fields[2]                  <-- ignore

                    if not fields[3].upper()  in ('DNA', 'RNA'):
                        raise Exception("'DNA' or 'RNA' expected, got '%s' instead.  Line:\n%s" % (fields[3],value))
                    sequence_record.molecule_type = fields[3]
                    fields=fields[4:]

                    # Last item on the line is the date:
                    try:
                        sequence_record.date = datetime.datetime.strptime(fields.pop(), '%d-%b-%Y')
                    except:
                        sequence_record.date = 'UNKNOWN'

                    # Division: optional, comes just before date
                    if len(fields) > 0 :
                       sequence_record.division = fields.pop()

                    """sequence_record.locus = fields[0]
                    sequence_type = fields[2]

                    last_field = fields[-1].split(' ')
                    if len(last_field) == 3:
                        sequence_record.molecule_type = last_field[0]

                    # Division

                    if len(last_field) > 1:
                        sequence_record.division = last_field[-2]

                    for division in ('PRI', 'ROD', 'MAM', 'VRT', 'INV', 'PLN',
                                     'BCT', 'VRL', 'PHG', 'SYN', 'UNA', 'EST',
                                     'PAT', 'STS', 'GSS', 'HTG', 'HTC', 'ENV'):
                        if division in fields:
                            sequence_record.division = division

                    sequence_record.date = datetime.datetime.strptime(last_field[-1],
                                                                      '%d-%b-%Y')

                    import pdb
                    pdb.set_trace()"""

                elif attribute == 'DEFINITION':
                    #print value
                    sequence_record.definition = value

                elif attribute == 'ACCESSION':
                    sequence_record.accession = value

                elif attribute == 'VERSION':
                    sequence_record.version = value

                elif attribute == 'KEYWORDS':
                    sequence_record.keywords = [keyword for keyword in
                                                re.split(r';\s*',
                                                         value.rstrip('.'))
                                                if keyword]

                elif attribute == 'SOURCE':
                    sequence_record.source = value

                    while lines and self.SUBATTRIBUTE_REGEXP.match(lines[0]):
                        line = lines.pop(0)

                        subattribute = line[2:11].strip()
                        subvalue = [line[12:]]

                        while lines and lines[0].startswith(' ' * 11):
                            subvalue.append(lines.pop(0)[12:])

                        if subattribute == 'ORGANISM':
                            hierarchy = re.split(r';\s*',
                                                 ''.join(subvalue[1:]).rstrip('.'))
                            sequence_record.source_organisms.append({'organism':
                                                                     subvalue[0],
                                                                     'taxonomy':
                                                                     hierarchy})

                elif attribute == 'REFERENCE':
                    reference_id = value
                    authors = []
                    consortium = None
                    title = None
                    journal = None
                    pubmed = None

                    while lines and self.SUBATTRIBUTE_REGEXP.match(lines[0]):
                        line = lines.pop(0)

                        subattribute = line[2:11].strip()
                        subvalue = [line[12:]]

                        while lines and lines[0].startswith(' ' * 11):
                            subvalue.append(lines.pop(0)[12:])

                        if subattribute == 'AUTHORS':
                            authors = re.split(r',\s+',
                                               ' '.join(subvalue).replace(' and ', ', '))

                        elif subattribute == 'CONSRTM':
                            consortium = ' '.join(subvalue)

                        elif subattribute == 'TITLE':
                            title = ' '.join(subvalue)

                        elif subattribute == 'JOURNAL':
                            journal = ' '.join(subvalue)

                        elif subattribute == 'PUBMED':
                            pubmed = ' '.join(subvalue)

                    sequence_record.references.append({'id': reference_id,
                                                       'authors': authors,
                                                       'consortium':
                                                       consortium,
                                                       'title': title,
                                                       'journal': journal,
                                                       'pubmed': pubmed})

                elif attribute == 'COMMENT':
                    comment = [value]
                    while lines and lines[0].startswith(' ' * 11):
                        comment.append(lines.pop(0)[12:])

                    sequence_record.comments.append('\n'.join(comment))

            # Discard empty lines
            if not lines[0].strip():
                del lines[0]

            if not lines:
                raise IndexError('Unexpected end of file encountered')
            if not lines[0].startswith('FEATURES'):
                raise ValueError('Unexpected token encountered while parsing features: %s'
                                 % repr(lines[0]))
            
            # Parse the features out
            del lines[0]
            while (lines and not lines[0].startswith('ORIGIN') and not
                   lines[0].startswith('BASE COUNT')):

                # Ignore empty lines

                if not lines[0].strip():
                    del lines[0]
                    continue

		# Match source, gene, CDS or anything starting with 5 spaces and a letter in column 6
                if self.FEATURE_REGEXP.match(lines[0]):
                    current_line = lines.pop(0).strip()

                    # Attempt to glob more lines, in case the location field
                    # extends beyond one line (rare)
                    while (lines and not lines[0].lstrip().startswith('/') and
                           lines[0].startswith(' ' * 7)):
                        current_line += lines.pop(0).strip()

                    # Parse out the feature type

                    feature_type, location_string = re.split(r'\s+',
                                                             current_line, 1)

                    # Parse out the location, for some reason

                    # The fact that there are more than one entries in location_strings[] implies 'join'.
                    if 'join(' in location_string:
                        location_string = location_string.replace('join(', '').rstrip(')')

                    location_strings = location_string.split(',')

                    locations = []
                    for location_string in location_strings:
                        start = None
                        end = None

                        if 'complement(' in location_string:
                            location_string = location_string.replace('complement(',
                                                                      '').rstrip(')')
                            try:
                                end, start = location_string.split('..')
                            except ValueError:
                                if re.match(r'^\d+$', location_string):
                                    start, end = location_string, location_string
                                else:
                                    raise ValueError('Could not parse the location %s from line:\n%s'
                                                     % (repr(location_string),
                                                        repr(current_line)))
                        else:
                            try:
                                start, end = location_string.split('..')
                            except ValueError:
                                if re.match(r'^\d+$', location_string):
                                    start, end = location_string, location_string
                                else:
                                    raise ValueError('Could not parse the location %s from line:\n%s'
                                                     % (repr(location_string),
                                                        repr(current_line)))

                        start_matches = self.FUZZY_LOCATION_REGEXP.match(start)
                        end_matches = self.FUZZY_LOCATION_REGEXP.match(end)

                        start_coordinate = (location.FuzzyCoordinate(start_matches.group(0),
                                                                     start_matches.group(1))
                                            if start_matches
                                            else location.Coordinate(start))
                        end_coordinate = (location.FuzzyCoordinate(end_matches.group(0),
                                                                     end_matches.group(1))
                                            if end_matches
                                            else location.Coordinate(end))
                        
                        locations.append(location.Coordinates(start_coordinate,
                                                              end_coordinate))

                    feature_location = location.Location(locations)
                    #feature_location = location_string

                    # Pull qualifiers out
                    qualifiers = {}
                    while lines and self.QUALIFIER_REGEXP.match(lines[0]):
                        parts = lines.pop(0).lstrip(' /').split('=', 1)

                        qualifier = parts[0]
                        value = [parts[1]] if len(parts) > 1 else []

                        while (lines and lines[0].startswith(' ' * 21)
                               #and not lines[0].strip().startswith('/')):
                               and not re.search(r'^/[A-Za-z_]+=', lines[0].strip())):
                            value.append(lines.pop(0).lstrip())
                        
                        # Special thing: if we're dealing with a translation,
                        # we will join with the empty string; otherwise, join
                        # with a space

                        if qualifier == 'translation':
                            value = ''.join(value)

                        # If we're dealing with a note field, we'll also have
                        # to do something special

                        elif qualifier == 'note':

                            individual_notes = []
                            current_note = []
                            while value:
                                current_line = value.pop(0)
                                if (re.match(r'^[A-Za-z]+:', current_line) and
                                    current_note):
                                    individual_notes.append(' '.join(current_note))
                                    current_note = []
                                current_note.append(current_line)
                                if not value:
                                    individual_notes.append(' '.join(current_note))
                            value = '\n'.join(individual_notes)
                                    
                        else:
                            value = ' '.join(value)

                        value = re.sub(r'\s{2,}', '\n', value.strip('"'))

                        # If no value was found, set it to None
                        if not value:
                            value = None

                        if qualifier not in qualifiers:
                            qualifiers[qualifier] = []
                        qualifiers[qualifier].append(value)

                    # Add a locus_tag to gene and CDS features that lack one:
                    if 'locus_tag' not in qualifiers and feature_type in ('gene', 'CDS'):
                        qualifiers['locus_tag'] = []
                        qualifiers['locus_tag'].append('_'.join((sequence_record.locus,'%04d' % ltnum)))   # Add a locus_tag
                        # Only increment the locus tag number for CDS entries.  Assumes that CDS follows corresponding gene feature.
                        if feature_type == 'CDS':
                            ltnum += 1

                    # Add this feature to the list of features:
                    sequence_record.features.append(Feature(feature_type,
                                                            feature_location,
                                                            qualifiers=qualifiers))

                else:
                    raise ValueError('Unexpected token encountered while parsing features: %s'
                                     % repr(lines[0]))


	    # At this point, only the nucleotide sequence remains:
            if lines and lines[0].startswith('BASE COUNT'):
                del lines[0]

            del lines[0]  # ORIGIN line deleted

            sequence=[]
            while lines:
                line = lines.pop(0)
                if not line.startswith('//'):    # If not the end
                    sequence.extend(line.strip().upper().split(' ')[1:])  # Add more sequence bits to the list ([0] = number)

            sequence_record.sequence=''.join(sequence)  # Make a big string out of it

            return sequence_record
        raise StopIteration()



class GenBankRecordSerializer(SequenceRecordSerializer):
    """Writes a record to the GenBank format."""
    ATTRIBUTE_WRAPPER = textwrap.TextWrapper(width=79,
                                             subsequent_indent=' ' * 12)
    QUALIFIER_WRAPPER = textwrap.TextWrapper(width=79,
                                             initial_indent=' ' * 21,
                                             subsequent_indent=' ' * 21)

    def __serialize_attribute(self, attribute, value):
        return '\n'.join(self.ATTRIBUTE_WRAPPER.wrap('%-11s %s' % (attribute,
                                                                   value)))

    def __serialize_feature(self, feature):
        return '     %-15s %s' % (feature.type, feature.location)

    def __serialize_qualifier(self, qualifier, value):
        if value:
            output_lines = []
            lines = value.split('\n')
            if len(lines) == 1:
                output_lines.append('\n'.join(self.QUALIFIER_WRAPPER.wrap('/%s="%s"' %
                                                                          (qualifier,
                                                                           lines.pop(0)))))
            else:
                output_lines.append('\n'.join(self.QUALIFIER_WRAPPER.wrap('/%s="%s' %
                                                                          (qualifier,
                                                                           lines.pop(0)))))
                while len(lines):
                    line = '\n'.join(self.QUALIFIER_WRAPPER.wrap(lines.pop(0)))
                    if lines:
                        output_lines.append(line)
                    else:
                        output_lines.append(line + '"')
            return '\n'.join(output_lines)
        else:
            return '\n'.join(self.QUALIFIER_WRAPPER.wrap('/%s' % qualifier))

    def serialize(self, record):
        lines = []

        # Locus
        #print >>sys.stderr,"record.date=%s" % record.date;
        lines.append('LOCUS       %-20s %7d bp    %-3s  %-3s %11s' %
                     (record.locus, 
                      len(record.sequence) if record.sequence else 0,
                      record.molecule_type if record.molecule_type else 'DNA',
                      record.division if record.division else '',
                      record.date.strftime('%d-%b-%Y').upper()))
        #print >>sys.stderr,"lines[-1]: %s" % lines[-1]
        #sys.exit(17);

        # Definition
        lines.append(self.__serialize_attribute('DEFINITION',
                                                record.definition))

        # Accession
        lines.append(self.__serialize_attribute('ACCESSION', record.accession))

        # Version
        lines.append(self.__serialize_attribute('VERSION', record.version))

        # Keywords
        lines.append(self.__serialize_attribute('KEYWORDS',
                                                '; '.join(record.keywords) +
                                                '.'))

        # Source
        lines.append(self.__serialize_attribute('SOURCE', record.source))

        # Source organisms
        for source_organism in record.source_organisms:
            lines.append(self.__serialize_attribute('  ORGANISM', source_organism['organism']))
            lines.append(self.__serialize_attribute('',
                                                    '; '.join(source_organism['taxonomy'])
                                                    + '.'))

        # References
        for reference in record.references:
            lines.append(self.__serialize_attribute('REFERENCE',
                                                    reference['id']))
            if reference['authors']:
                if len(reference['authors']) > 1:
                    lines.append(self.__serialize_attribute('  AUTHORS', '%s and %s'
                                                            % (', '.join(reference['authors'][0:-1]),
                                                               reference['authors'][-1])))
                else:
                    lines.append(self.__serialize_attribute('  AUTHORS', reference['authors'][0]))

            if reference['consortium']:
                lines.append(self.__serialize_attribute('  CONSRTM',
                                                        reference['consortium']))

            if reference['title']:
                lines.append(self.__serialize_attribute('  TITLE',
                                                        reference['title']))

            if reference['journal']:
                lines.append(self.__serialize_attribute('  JOURNAL',
                                                        reference['journal']))

            if reference['pubmed']:
                lines.append(self.__serialize_attribute('  PUBMED',
                                                        reference['pubmed']))

        # Comment

        for comment in record.comments:
            comment_lines = comment.split('\n')
            lines.append(self.__serialize_attribute('COMMENT',
                                                    comment_lines.pop(0)))
            for comment_line in comment_lines:
                lines.append(self.__serialize_attribute('', comment_line))

        # Features

        lines.append('FEATURES             Location/Qualifiers')

        for feature in record.features:
            lines.append(self.__serialize_feature(feature))
            for qualifier, values in feature.qualifiers.items():

                # db_xref: sort the values

                if qualifier == 'db_xref':
                    values = sorted(values)

                for value in values:
                    lines.append(self.__serialize_qualifier(qualifier, value))

        # Sequence

        lines.append('ORIGIN')

        if record.sequence:
            i = 1
            parts = []
            part = ''
            for n in iter(record.sequence):
                part += n
                if len(part) == 10:
                    parts.append(part)
                    part = ''
                if len(parts) == 6:
                    lines.append('%9d %s' % (i, ' '.join(parts)))
                    parts = []
                    i += 60
            if part:
                parts.append(part)
            if parts:
                lines.append('%9d %s' % (i, ' '.join(parts)))

        lines.append('//')

        return '\n'.join(lines)
