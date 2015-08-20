# py-vcf-parser

Simple [VCF (variant call format)](https://github.com/samtools/hts-specs) parser for python.


## Getting Started

```
$ pip install py-vcf-parser
```


## Usage example

```python
>>> from vcf.reader import VCFReader
>>> with open('test/test.vcf41.vcf', 'r') as fin:
...     for record in VCFReader(fin):n
...         print record
{'ALT': ['A'],
 'CHROM': '20',
 'FILTER': 'PASS',
 'ID': 'rs6054257',
 'INFO': 'NS=3;DP=14;AF=0.5;DB;H2',
 'NA00001': {'DP': '1', 'GQ': '48', 'GT': '0|0', 'HQ': '51,51'},
 'NA00002': {'DP': '8', 'GQ': '48', 'GT': '1|0', 'HQ': '51,51'},
 'NA00003': {'DP': '5', 'GQ': '43', 'GT': '1/1', 'HQ': '.,.'},
 'QUAL': '29',
 'REF': 'G',
 'allele_count': 6,
 'allele_freq': {'A': 0.5, 'G': 0.5},
 'chrom': '20',
 'genotype': {'NA00001': ['G', 'G'],
              'NA00002': ['A', 'G'],
              'NA00003': ['A', 'A']},
 'info': {'AF': '0.5', 'DB': True, 'DP': '14', 'H2': True, 'NS': '3'},
 'pos': 14370,
 'rs': 6054257}
...
```


## Tests

```
$ pip install pytest
$ py.test
```

- Tested on python2.7 and pypy2 only. python3.x is not supported.
- Tested on VCF version 4.1.


## Notes

- `symbolic allele` and `breakend` are not supported.


## License

See `LICENSE.txt`


## Author

[@knmkr](https://github.com/knmkr)
