# -*- coding: utf-8 -*-

import sys
import csv
csv.field_size_limit(sys.maxsize)  # FIXME
import re


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
            data['genotype'] = {}

            data['QUAL'] = str(record['QUAL'])
            data['FILTER'] = str(record['FILTER'])

            data['INFO'] = str(record['INFO'])
            data['info'] = _info(record['INFO'])

            # FORMAT
            format_keys = record['FORMAT'].split(':')
            for sample in self.sample_names:
                data[sample] = dict(zip(format_keys, record[sample].split(':')))
                data['genotype'][sample] = _GT2genotype(data['REF'],
                                                        data['ALT'],
                                                        data[sample]['GT'])

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
    >>> _alt('G,T,A')
    ['G', 'T', 'A']
    """
    alt = text.split(',')
    if len(alt) != 1:
        print >>sys.stderr, "[WARN] multiple alts:", alt

    return alt

def _rsid(text):
    """
    >>> _rsid('rs100')
    100
    """

    # TODO: simply split by `;` maybe buggy...
    rs_raw = text.split(';')  # rs100;rs123
    rs_first = rs_raw[0]
    if len(rs_raw) != 1:
        print >>sys.stderr, "[WARN] multiple rsids:", rs_raw

    rs_regex = re.compile('rs(\d+)')
    rs_match = rs_regex.match(rs_first)

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

    >>> info = _info('NS=3;DP=14;AF=0.5;DB;H2')
    >>> info['H2']
    True
    >>> info['NS']
    '3'
    >>> info['DB']
    True
    >>> info['DP']
    '14'
    >>> info['AF']
    '0.5'
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
    'GG'
    >>> GT = '0'  # 1 allele (chrX, etc.)
    >>> _GT2genotype(REF, ALT, GT)
    'G'
    """

    # / : genotype unphased
    # | : genotype phased

    gt = re.split('/|\|', GT)

    # 0 : reference allele (what is in the REF field)
    # 1 : first allele listed in ALT
    # 2 : second allele list in ALT
    # and so on.

    bases = [REF] + ALT

    if len(gt) == 2:
        genotype = bases[int(gt[0])] + bases[int(gt[1])]
    elif len(gt) == 1:
        genotype = bases[int(gt[0])]

    return genotype
