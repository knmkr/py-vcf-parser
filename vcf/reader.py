# -*- coding: utf-8 -*-

import sys
import csv
csv.field_size_limit(sys.maxsize)  # FIXME
import re
from collections import Counter
from pprint import pprint
from decimal import *

class VCFFilterNotPass(Exception):
    pass

class VCFReader(object):
    def __init__(self, fin, filters={}, decimal_prec=4):
        self.fin = fin
        self.delimiter = '\t'
        self.headerlines = []
        self.fieldnames = []

        # Header lines
        for line in self.fin:
            line = line.rstrip()

            if line.startswith('##'):
                self.headerlines.append(line)
                continue
            elif line.startswith('#CHROM'):
                self.fieldnames = line.replace('#CHROM', 'CHROM').split(self.delimiter)
                break
            else:
                raise csv.Error, 'Invalid header lines. `#CHROM ...` does not exists.'

        if self.fieldnames[0:9] != ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']:
            raise csv.Error, 'Invalid header lines. Probably delimiter is not tab.'

        self.sample_names = self.fieldnames[9:]
        self.filters = filters
        self.decimal_prec = decimal_prec

    def __iter__(self):
        data = {}
        for record in csv.DictReader(self.fin,
                                     fieldnames=self.fieldnames,
                                     delimiter=self.delimiter):
            try:
                data['CHROM'] = str(record['CHROM'])
                data['chrom'] = _chrom(record['CHROM'])

                data['pos'] = int(record['POS'])

                data['ID'] = str(record['ID'])
                data['rs'] = _rsid(record['ID'])

                data['REF'] = str(record['REF'])
                data['ALT'] = _alt(record['ALT'])

                data['QUAL'] = str(record['QUAL'])
                data['FILTER'] = str(record['FILTER'])

                data['INFO'] = str(record['INFO'])
                data['info'] = _info(record['INFO'])

                data['genotype'] = {}
                for sample in self.sample_names:
                    data[sample] = dict(zip(record['FORMAT'].split(':'), _sample(record[sample]).split(':')))
                    data['genotype'][sample] = _GT2genotype(data['REF'],
                                                            data['ALT'],
                                                            data[sample]['GT'])

                for key,condition in self.filters.items():
                    if key == 'genotype':
                        data[key] = dict(filter(condition, data[key].iteritems()))
                    elif key == 'rs':
                        if not condition(data[key]):
                            raise VCFFilterNotPass

                counter = count_allele(data['genotype'])
                data['allele_count'] = sum(counter.values())

                getcontext().prec = self.decimal_prec
                data['allele_freq'] = {k:allele_freq(cnt, data['allele_count']) for k,cnt in counter.items()}

                if data['allele_count'] == 0:
                    print >>sys.stderr, '[WARN] allele_count is 0:', data['ID']

            except VCFFilterNotPass:
                continue
            else:
                yield data


def _chrom(text):
    if text.startswith('chr'):
        text = text.replace('chr', '')

    if text in [str(i+1) for i in range(22)] + ['X', 'Y', 'MT', 'M']:
        return text
    else:
        return None

def _alt(text):
    """
    >>> _alt('.')
    ['.']
    >>> _alt('G')
    ['G']
    >>> _alt('G,T')
    ['G', 'T']
    >>> _alt('G,T,A')
    ['G', 'T', 'A']
    >>> _alt('G,T,A,C')
    ['G', 'T', 'A', 'C']
    >>> _alt('<INS:ME:ALU>')
    ['<INS:ME:ALU>']
    >>> _alt('FOO')
    Traceback (most recent call last):
        ...
    Error: Invalid ALT field: FOO
    """
    regex = re.compile('([ACGTN,.]+|<.+>)')
    if not regex.match(text):
        raise csv.Error, 'Invalid ALT field: {ALT}'.format(ALT=text)

    alt = text.split(',')

    return alt

def _rsid(text):
    """
    >>> _rsid('rs100')
    100
    >>> type(_rsid('rs100')) == int
    True
    >>> _rsid('.')
    """

    # TODO: handle multiple rs like rs100;rs123
    regex = re.compile('rs(\d+)')
    rs_match = regex.match(text)

    if rs_match:
        rsid = rs_match.group(1)
    else:
        return None

    try:
        return int(rsid)
    except ValueError:
        return None

def _info(text):
    """parse INFO field.

    >>> pprint(_info('NS=3;DP=14;AF=0.5;DB;H2'))
    {'AF': '0.5', 'DB': True, 'DP': '14', 'H2': True, 'NS': '3'}
    """

    info = dict()
    for record in text.split(";"):
        try:
            k,v = record.split("=")
        except ValueError:
            k,v = (record, True)
        info.update({k:v})
    return info

def _sample(text):
    if text == None:
        return ''
    else:
        return text

def _GT2genotype(REF, ALT, GT):
    """parse GT (GenoType) in Genotype fields.

    >>> _GT2genotype('G', ['A'], '0|0')
    ['G', 'G']
    >>> _GT2genotype('G', ['A'], '0/0')  # unphased genotype
    ['G', 'G']
    >>> _GT2genotype('G', ['A'], '0')    # 1 allele (chrX, etc.)
    ['G']
    >>> _GT2genotype('G', ['A'], './.')  # no calls
    []
    >>> _GT2genotype('G', ['A'], '')     # blank
    []
    >>> _GT2genotype('G', ['A','T'], '1|2')
    ['A', 'T']
    >>> _GT2genotype('G', ['A','T','C'], '3|3')
    ['C', 'C']
    """

    if GT == '':
        return []

    # / : genotype unphased
    # | : genotype phased

    gt = re.split('/|\|', GT)

    # 0 : reference allele (what is in the REF field)
    # 1 : first allele listed in ALT
    # 2 : second allele list in ALT
    # and so on.

    alleles = [REF] + ALT
    genotype = [alleles[int(idx)] for idx in gt if idx != '.']

    return genotype

def count_allele(genotype):
    """
    >>> cnt = count_allele({'N1': ['A', 'A'], 'N2': ['A', 'T'], 'N3': ['T', 'T']})
    >>> cnt['A']
    3
    >>> cnt['T']
    3
    >>> cnt['G']
    0
    >>> cnt = count_allele({'N1': ['G', 'G'], 'N2': ['G', 'G'], 'N3': ['G', 'G']})
    >>> cnt['G']
    6
    """

    counter = Counter()
    for sample, alleles in genotype.items():
        for allele in alleles:
            counter[allele] += 1
    return counter

def allele_freq(count, total):
    """
    >>> getcontext().prec = 3
    >>> allele_freq(1,3)
    Decimal('0.333')
    >>> allele_freq(2675,10000)
    Decimal('0.268')
    >>> getcontext().prec = 4
    >>> allele_freq(2675,10000)
    Decimal('0.2675')
    >>> allele_freq(2675,100000)
    Decimal('0.02675')
    """

    return Decimal(count) / Decimal(total)
