# -*- coding: utf-8 -*-

import sys
import csv
csv.field_size_limit(sys.maxsize)  # FIXME


class VCFWriter(object):
    def __init__(self, fout, sample_ids, header_records):
        self.fout = fout
        self.sample_ids = sample_ids
        self.header_records = header_records

        fieldnames = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + self.sample_ids
        self.writer = csv.DictWriter(self.fout, fieldnames=fieldnames, delimiter='\t')

    def writeheaderlines(self):
        for record in self.header_records:
            self.fout.write('##{}={}\n'.format(record[0], record[1]))
        self.writer.writeheader()

    def writerow(self, record):
        info = []
        for x in record['INFO']:
            if x[1] == True:
                info.append(x[0])
            else:
                info.append('='.join(x))

        row =  {'#CHROM':          record['#CHROM'],
                'POS':             record['POS'],
                'ID':              record['ID'],
                'REF':             record['REF'],
                'ALT':    ','.join(record['ALT']),
                'QUAL':            record['QUAL'],
                'FILTER':          record['FILTER'],
                'INFO':   ';'.join(info),
                'FORMAT': ':'.join(record['FORMAT'])}

        for sample_id in self.sample_ids:
            row[sample_id] = ':'.join([x[1] for x in record[sample_id]])

        self.writer.writerow(row)
