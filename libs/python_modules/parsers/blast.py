"""This module defines classes for parsing BLAST output."""

import multiprocessing
import time
import sys
import os
import re
import math

def parse_entry(entry, options):
    # Pull out details about the entry

    query = None
    lines = entry.split("\n")
    results = []
    while lines:

        line = lines.pop(0)

        # Set the name of our query if applicable
        
        if line.startswith('Query='):
            query = line[7:].strip()

        # Start pulling if our line starts with '>'

        if line.startswith('>'):
            lines.insert(0, line)
            result = dict(query=query, kegg=None, ec_numbers=[],
                          match_length=0, source_length=0, hit_coverage=0.0,
                          score=0.0, expect=1e1000, identities=0.0,
                          positives=0.0, gaps=100.0)

            # Parse the product

            product = []
            while lines and 'Length' not in lines[0]:
                product.append(lines.pop(0).strip())
            product = ' '.join(product).lstrip('>')

            # If KEGG mode is enabled, pull out the identifier

            if options['kegg_mode']:
                kegg_results = re.search(r'^([a-z]{3}:[A-Z0-9_]+)', product)
                if kegg_results:
                    options['kegg'] = kegg_results.group(1)
                    product = re.sub(r'^([a-z]{3}:[A-Z0-9_]+)', '', product)
            result['product'] = product.strip()

            # If there are any EC numbers, pull them out!

            ec_results = re.findall(r'\[EC:(.+?)\]', product)
            result['ec_numbers'] = []
            for ec_result in ec_results:
                result['ec_numbers'].extend(re.split(r'\s+', ec_result))

            # Pull out the length of the source

            if lines and 'Length' in lines[0]:
                result['source_length'] = int(lines.pop(0).strip().split(' = ')[1])

            while lines and 'Score' not in lines[0]:
                lines.pop(0)

            # Pull out statistics

            if lines and 'Score' in lines[0]:
                subline = lines.pop(0)

                score_results = re.search(r'Score =\s+(\d+(\.\d+)?)', subline)
                if score_results:
                    result['score'] = float(score_results.group(1))

                    # Skip if the score is too low
                    if result['score'] < options['minimum_score']:
                        continue

                expect_results = re.search(r'Expect =\s+(\S+)', subline)
                if expect_results:
                    expect_value = expect_results.group(1)
                    if expect_value.startswith('e'):
                        expect_value = '1%s' % expect_value
                    result['expect'] = float(expect_value)

                    # Skip if the expect value is too high
                    if result['expect'] > options['maximum_expect_value']:
                        continue

            if lines and 'Identities' in lines[0]:
                subline = lines.pop(0)

                identities_result = re.search(r'Identities =\s+(\d+)/(\d+)', subline)
                if identities_result:
                    result['identities'] = (float(identities_result.group(1)) /
                                            float(identities_result.group(2)) *
                                            100)
                    result['identities_num'] = int(identities_result.group(1))
                    
                    # Skip if the % identity is too low
                    if result['identities'] < options['minimum_identity']:
                        continue

                    result['match_length'] = int(identities_result.group(2))

                    # Skip if the match length is too short or too long
                    if (result['match_length'] < options['minimum_length'] or
                        result['match_length'] > options['maximum_length']):
                        continue

                    result['hit_coverage'] = (result['match_length'] /
                                              float(result['source_length']) *
                                              100)

                    # Skip if the coverage is too low
                    if result['hit_coverage'] < options['minimum_hit_coverage']:
                        continue

                positives_result = re.search(r'Positives =\s+(\d+)/(\d+)', subline)
                if positives_result:
                    result['positives'] = (float(positives_result.group(1)) /
                                           float(positives_result.group(2)) *
                                           100)
                    result['positives_num'] = int(positives_result.group(1))

                    # Skip if the % positives is too low
                    if result['positives'] < options['minimum_positives']:
                        continue

                gaps_result = re.search(r'Gaps =\s+(\d+)/(\d+)', subline)
                if gaps_result:
                    result['gaps'] = (float(gaps_result.group(1)) /
                                      float(gaps_result.group(2)) * 100)
                    result['gaps_num'] = int(positives_result.group(1))

                    # Skip if the % gaps is too high
                    if result['gaps'] > options['maximum_gaps']:
                        continue

            while lines and not lines[0].startswith('>'):
                lines.pop(0)

            # Calculate the heuristic
            # The heuristic is: bit score * hit coverage / 100 * identities /
            # 100 * positives / 100 * (1 - gaps / 100) * max(100,
            # -log10(expect) / 100)
            # TODO: Tweak this heuristic
            heuristic = (result['score'] * result['hit_coverage'] / 100 *
                         result['identities'] / 100 * result['positives'] / 100 *
                         (1 - (result['gaps'] / 100)) * max(100, -1 *
                                                            math.log10(result['expect']))
                         / 100)
            result['heuristic'] = heuristic

            results.append(result)

    # Pick the top n results
    results.sort(cmp=lambda x, y: cmp(y['heuristic'], x['heuristic']))
    n = 0
    ret = []
    while n < int(options['number_of_results']) and results:
        result = results.pop(0)
        ret.append((query, result['product'], result['score'],
                    result['expect'], result['match_length'],
                    result['source_length'], result['hit_coverage'],
                    result['identities_num'], result['identities'],
                    result['positives_num'], result['positives'],
                    result['gaps_num'], result['gaps'], result['ec_numbers'],
                    result['heuristic']))
        n += 1

    #sys.stderr.write('.')

    return ret



class BlastParser(object):
    """Parses BLAST output."""

    def __init__(self, handler, number_of_results=1, minimum_score=0.0,
                 minimum_hit_coverage=0.0, maximum_expect_value=1e-6,
                 minimum_length=0, maximum_length=1e1000, minimum_identity=0.0,
                 maximum_identity=100.0, minimum_positives=0.0,
                 maximum_gaps=100.0, bsr_file=None, minimum_bsr=0.0,
                 kegg_mode=False):
        self.handler = handler
        self.options = {
            'number_of_results': int(number_of_results),
            'minimum_score': float(minimum_score),
            'minimum_hit_coverage': float(minimum_hit_coverage),
            'maximum_expect_value': float(maximum_expect_value),
            'minimum_length': int(minimum_length),
            'maximum_length': float(maximum_length),
            'minimum_identity': float(minimum_identity),
            'maximum_identity': float(maximum_identity),
            'minimum_positives': float(minimum_positives),
            'maximum_gaps': float(maximum_gaps),
            'bsr_file': bsr_file,
            'minimum_bsr': float(minimum_bsr),
            'kegg_mode': bool(kegg_mode)
        }

    def __get_entry_line(self):
        line = self.handler.readline()
        if not line:
            return None
        if line.startswith('BLAST'):
            self.handler.seek(self.handler.tell() - len(line))
            return None
        return line

    def __get_results(self):
        results = []

        # Hello, pool

        pool = multiprocessing.Pool(processes=1)

        def worker_callback(result):
            results.extend(result)

        # We'll have to manually manipulate the file pointer here...

        while True:
            line = self.handler.readline()
            if line == '':
                break

            # Start pulling out an entry

            if line.startswith('BLAST'):
                entry = [line]
                entry_line = self.__get_entry_line()
                while entry_line is not None:
                    entry.append(entry_line)
                    entry_line = self.__get_entry_line()

                # Perform a dispatch to parse the entry

                pool.apply_async(parse_entry, args=[''.join(entry),
                                                    self.options],
                                 callback=worker_callback)

        pool.close()
        pool.join()

        return results

    results = property(__get_results)
