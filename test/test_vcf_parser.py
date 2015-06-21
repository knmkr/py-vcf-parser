# -*- coding: utf-8 -*-

import os
import csv
import unittest

import pytest

import vcf


class SimpleTest(unittest.TestCase):
    def setUp(self):
        self.basedir = os.path.dirname(os.path.abspath(__file__))

    def test_success(self):
        with open(os.path.join(self.basedir, 'test.vcf41.vcf'), 'r') as fin:
            reader  = vcf.DictReader(fin)

            assert reader.sample_names == ['NA00001', 'NA00002', 'NA00003']

            for record in reader:
                assert record['chrom'] == '20'
                assert record['pos'] == 14370
                assert record['ID'] == 'rs6054257'
                assert record['rs'] == 6054257
                assert record['REF'] == 'G'
                assert record['ALT'] == ['A']

                assert record['INFO'] == 'NS=3;DP=14;AF=0.5;DB;H2'
                assert record['info'] == {'H2': True, 'NS': '3', 'DB': True, 'DP': '14', 'AF': '0.5'}

                assert record['NA00001'] == {'GT': '0|0',
                                             'GQ': '48',
                                             'DP': '1',
                                             'HQ': '51,51'}
                assert record['NA00002'] == {'GT': '1|0',
                                             'GQ': '48',
                                             'DP': '8',
                                             'HQ': '51,51'}
                assert record['NA00003'] == {'GT': '1/1',
                                            'GQ': '43',
                                             'DP': '5',
                                             'HQ': '.,.'}

                assert record['genotypes'] == {'NA00001': 'GG',
                                               'NA00002': 'AG',
                                               'NA00003': 'AA'}

                assert record['samples'] == ['NA00001', 'NA00002', 'NA00003']

                break

    def test_header_without_chrom_should_fail_parse(self):
        with pytest.raises(csv.Error) as e:
            reader = vcf.DictReader(open(os.path.join(self.basedir, 'test.vcf41.invalid-header.header-without-chrom.vcf'), 'r'))

        assert 'Invalid header lines. `#CHROM ...` does not exists.' in str(e)


    def test_without_header_lines_should_fail_parse(self):
        with pytest.raises(csv.Error) as e:
            reader = vcf.DictReader(open(os.path.join(self.basedir, 'test.vcf41.invalid-header.without-header-lines.vcf'), 'r'))

        assert 'Invalid header lines. `#CHROM ...` does not exists.' in str(e)

    def test_delimiter_is_not_tab_should_fail_parse(self):
        with pytest.raises(csv.Error) as e:
            reader = vcf.DictReader(open(os.path.join(self.basedir, 'test.vcf41.invalid-header.delimiter-is-not-tab.vcf'), 'r'))

        assert 'Invalid header lines. Probably delimiter is not tab.' in str(e)

    def test_blank_file_should_fail_parse(self):
        with pytest.raises(csv.Error) as e:
            reader = vcf.DictReader(open(os.path.join(self.basedir, 'test.vcf.invalid.blank-file.vcf'), 'r'))

        assert 'Invalid header lines. Probably delimiter is not tab.' in str(e)


if __name__ == '__main__':
    unittest.main()
