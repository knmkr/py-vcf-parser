# -*- coding: utf-8 -*-

import sys
import csv
csv.field_size_limit(sys.maxsize)  # FIXME
import re
from collections import Counter

class DictReader(object):
    def __init__(self, fin):
        self.fin = fin
        self.delimiter = '\t'
        self.headerlines = []
        self.fieldnames = []

        if type(fin) != file:
            raise csv.Error, 'type(fin) is not file.'

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
        else:
            self.sample_names = self.fieldnames[9:]

    def __iter__(self):
        data = {}
        for record in csv.DictReader(self.fin,
                                     fieldnames=self.fieldnames,
                                     delimiter=self.delimiter):
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
                data[sample] = dict(zip(record['FORMAT'].split(':'), record[sample].split(':')))
                data['genotype'][sample] = _GT2genotype(data['REF'],
                                                        data['ALT'],
                                                        data[sample]['GT'])

            counter = count_allele(data['genotype'])
            data['allele_count'] = sum(counter.values())
            data['allele_freq'] = {k:float(v)/data['allele_count'] for k,v in counter.items()}

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
    >>> _alt('')
    ['']
    >>> _alt('G')
    ['G']
    >>> _alt('G,T')
    ['G', 'T']
    >>> _alt('G,T,A')
    ['G', 'T', 'A']
    >>> _alt('G,T,A,C')
    ['G', 'T', 'A', 'C']
    """
    alt = text.split(',')

    return alt

def _rsid(text):
    """
    >>> _rsid('rs100')
    100
    >>> type(_rsid('rs100')) == int
    True
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

    >>> _info('NS=3;DP=14;AF=0.5;DB;H2')
    {'H2': True, 'NS': '3', 'DB': True, 'DP': '14', 'AF': '0.5'}
    """

    info = dict()
    for record in text.split(";"):
        try:
            k,v = record.split("=")
        except ValueError:
            k,v = (record, True)
        info.update({k:v})
    return info

def _GT2genotype(REF, ALT, GT):
    """parse GT (GenoType) in Genotype fields.

    >>> REF = 'G'
    >>> ALT = ['A']
    >>> GT = '0|0'
    >>> _GT2genotype(REF, ALT, GT)
    ['G', 'G']
    >>> GT = '0'  # 1 allele (chrX, etc.)
    >>> _GT2genotype(REF, ALT, GT)
    ['G']
    >>> REF = 'G'
    >>> ALT = ['A','T']
    >>> GT = '1|2'
    >>> _GT2genotype(REF, ALT, GT)
    ['A', 'T']
    """

    # / : genotype unphased
    # | : genotype phased

    gt = re.split('/|\|', GT)

    # 0 : reference allele (what is in the REF field)
    # 1 : first allele listed in ALT
    # 2 : second allele list in ALT
    # and so on.

    alleles = [REF] + ALT

    if len(gt) == 2:
        genotype = [alleles[int(gt[0])], alleles[int(gt[1])]]
    elif len(gt) == 1:
        genotype = [alleles[int(gt[0])]]
    else:
        raise csv.Error, 'Invalid GT (Genotype) field. len(gt) should be 1 or 2: {GT}'.format(GT=GT)

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
