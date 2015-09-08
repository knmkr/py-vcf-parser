# -*- coding: utf-8 -*-

import sys
import os
import csv
import unittest
import tempfile

import pytest

import vcf


class VCFWriterTest(unittest.TestCase):
    def setUp(self):
        self.basedir = os.path.dirname(os.path.abspath(__file__))

    def test_success(self):
        header_records = [('fileformat', 'VCFv4.1'),
                          ('reference', 'file:///seq/references/1000GenomesPilot-NCBI36.fasta'),
                          ('contig', '<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>'),
                          ('FORMAT', '<ID=GT,Number=1,Type=String,Description="Genotype">'),
                          ('FORMAT', '<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'),
                          ('FORMAT', '<ID=DP,Number=1,Type=Integer,Description="Read Depth">'),
                          ('FORMAT', '<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">')]

        sample_ids = ['NA00001', 'NA00002', 'NA00003']

        fout = tempfile.TemporaryFile(mode='w+t')
        writer = vcf.writer.VCFWriter(fout, sample_ids, header_records)
        writer.writeheaderlines()
        writer.writerow({'#CHROM': '20',
                         'POS': 14370,
                         'ID': 'rs6054257',
                         'REF': 'G',
                         'ALT': ['A'],
                         'QUAL': 29,
                         'FILTER': 'PASS',
                         'INFO': [('NS', '3'), ('DP', '14'), ('AF', '0.5'), ('DB', True), ('H2', True)],
                         'FORMAT': ['GT', 'GQ', 'DP', 'HQ'],
                         'NA00001': [('GT', '0|0'),
                                     ('GQ', '48'),
                                     ('DP', '1'),
                                     ('HQ', '51,51')],
                         'NA00002': [('GT', '1|0'),
                                     ('GQ', '48'),
                                     ('DP', '8'),
                                     ('HQ', '51,51')],
                         'NA00003': [('GT', '1/1'),
                                     ('GQ', '43'),
                                     ('DP', '5'),
                                     ('HQ', '.,.')]})

        fout.seek(0)
        lines = [line.strip() for line in fout.readlines()]
        assert lines == ['##fileformat=VCFv4.1',
                         '##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta',
                         '##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>',
                         '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                         '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
                         '##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">',
                         '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003',
                         '20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.']


if __name__ == '__main__':
    unittest.main()
