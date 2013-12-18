"""This module defines classes for working with KEGG ko files."""

import re

class KOEntry(object):
    """A ko entry."""
    def __init__(self, entry_id, name, definition, ec_numbers, classes,
                 database_links, genes):
        self.id = entry_id
        self.name = name
        self.definition = definition
        self.ec_numbers = ec_numbers
        self.classes = classes if classes else {}
        self.database_links = database_links if database_links else {}
        self.genes = genes if genes else {}

    def __repr__(self):
        return ('<KOEntry %s with %d classes, %d databases, and %d organisms>'
                % (self.id, len(self.classes), len(self.database_links),
                   len(self.genes)))



class KOParser(object):
    """A ko file parser."""
    def __init__(self, contents):
        self.contents = contents
        self.entries = []

    def __iter__(self):
        self.entries = [entry for entry in self.contents.split('///\n')
                        if entry]
        return self

    def next(self):
        if self.entries:
            entry = self.entries.pop(0)
            lines = [line for line in entry.split('\n') if line]
            
            # Parse out the features...
            features = {}

            while lines:
                first_line = lines.pop(0)
                feature_name, feature = first_line[0:12].strip(), [first_line]
                while lines and lines[0].startswith(' ' * 12):
                    feature.append(lines.pop(0))
                features[feature_name] = '\n'.join([l[12:] for l in feature])
            
            # Parse the entry feature out
            if 'ENTRY' in features:
                features['ENTRY'] = features['ENTRY'].split(' ', 1)[0]

            # Pull EC numbers out
            ec_numbers = []
            definition = None
            if 'DEFINITION' in features:
                matches = re.search(r'\s+\[EC:(.+)\]$', features['DEFINITION'])
                if matches:
                    ec_numbers = re.split(r'\s+', matches.group(1))
                    definition = features['DEFINITION'].replace(matches.group(1), '')
                else:
                    definition = features['DEFINITION']

            # Pull the classes out
            classes = []
            if 'CLASS' in features:
                class_lines = features['CLASS'].split('\n')
                while class_lines:
                    current_class = class_lines.pop(0)
                    while class_lines and not current_class.endswith(']'):
                        current_class += ' %s' % class_lines.pop(0)
                    matches = re.match(r'^(.+?)(\s+\[(.+)\])?$', current_class)
                    classes.append((matches.group(1), matches.group(3)))

            # Deal with database links
            database_links = {}
            if 'DBLINKS' in features:
                db_lines = features['DBLINKS'].split('\n')
                while db_lines:
                    db_line = db_lines.pop(0)
                    while db_lines and db_lines[0].startswith(' '):
                        db_line += ' %s' % db_lines.pop(0)
                    database, links = db_line.split(': ', 1)
                    links = links.split(' ')
                    database_links[database] = links

            # Deal with genes
            genes = {}
            if 'GENES' in features:
                gene_lines = features['GENES'].split('\n')
                while gene_lines:
                    gene_line = gene_lines.pop(0)
                    while gene_lines and gene_lines[0].startswith(' '):
                        gene_line += ' %s' % gene_lines.pop(0).strip()
                    organism, references = gene_line.split(': ', 1)
                    references = references.split(' ')
                    genes[organism] = references

            return KOEntry(features['ENTRY'], features['NAME'], definition,
                           ec_numbers, classes, database_links, genes)

        raise StopIteration()
