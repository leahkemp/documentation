# Compare dbSNP files/builds

Created: 2020/03/19 14:18:43
Last modified: 2020/03/26 13:36:06

## Table of contents

- [Compare dbSNP files/builds](#compare-dbsnp-filesbuilds)
  - [Table of contents](#table-of-contents)
  - [Findings](#findings)
  - [Comparison](#comparison)
    - [GRCh37](#grch37)
      - [Data sources](#data-sources)
    - [GRCh38](#grch38)
      - [Data sources](#data-sources-1)

This document compares the chromosome labelling of dbSNP databases of differing builds and downloaded from different sources. This can be useful when mitigating errors arising from chromosome labelling when building/running genomic pipelines that utilise these resources.

## Findings

- The chromosome labelling of the same dbSNP build is the same between GRCh37 and GRCh38
- The chromosome labelling differs between dbSNP databases hosted by NCBI and GATK, although NCBI provides a version of the dbSNP database that is compatible with GATK
- The chromosome labelling of the (current) newest dbSNP build (153) is different from older builds, instead using genbank accession numbers (see [here](https://www.ncbi.nlm.nih.gov/grc/human))
- The (current) newest dbSNP build (153) is only available from NCBI
- GATK currently only hosts older versions of dbSNP databases (build 138 for GRCh17 and build 146 for GRCh38)
- Some dbSNP databases provide the associated tabix index file (.tbi), however some do not and will need to be manually created with [tabix](http://www.htslib.org/doc/tabix.html)
- Some dbSNP databases are in the wrong compression format when downloaded and will need to be unzipped with [gunzip](https://linux.die.net/man/1/gunzip), and rezipped into bgzip format with [bgzip](http://www.htslib.org/doc/bgzip.html) before tabix can be used to create a tabix index file (.tbi)

Code used to return unique values of the chromosome column:

```bash
zgrep -v "#" file.gz | awk '{print $1}' | uniq
```

## Comparison

### GRCh37

| GATK (build 138)                   | NCBI (build 150) | NCBI (build 150) gatk version | NCBI (build151) | NCBI (build151) gatk version | NCBI (build 153) |
|------------------------------------|------------------|-------------------------------|-----------------|------------------------------|------------------|
| chrM                               | 1                | chr1                          | 1               | chr1                         | NC_000001.10     |
| chr1                               | 2                | chr2                          | 2               | chr2                         | NC_000002.11     |
| chr2                               | 3                | chr3                          | 3               | chr3                         | NC_000003.11     |
| chr3                               | 4                | chr4                          | 4               | chr4                         | NC_000004.11     |
| chr4                               | 5                | chr5                          | 5               | chr5                         | NC_000005.9      |
| chr5                               | 6                | chr6                          | 6               | chr6                         | NC_000006.11     |
| chr6                               | 7                | chr7                          | 7               | chr7                         | NC_000007.13     |
| chr7                               | 8                | chr8                          | 8               | chr8                         | NC_000008.10     |
| chr8                               | 9                | chr9                          | 9               | chr9                         | NC_000009.11     |
| chr9                               | 10               | chr10                         | 10              | chr10                        | NC_000010.10     |
| chr10                              | 11               | chr11                         | 11              | chr11                        | NC_000011.9      |
| chr11                              | 12               | chr12                         | 12              | chr12                        | NC_000012.11     |
| chr12                              | 13               | chr13                         | 13              | chr13                        | NC_000013.10     |
| chr13                              | 14               | chr14                         | 14              | chr14                        | NC_000014.8      |
| chr14                              | 15               | chr15                         | 15              | chr15                        | NC_000015.9      |
| chr15                              | 16               | chr16                         | 16              | chr16                        | NC_000016.9      |
| chr16                              | 17               | chr17                         | 17              | chr17                        | NC_000017.10     |
| chr17                              | 18               | chr18                         | 18              | chr18                        | NC_000018.9      |
| chr18                              | 19               | chr19                         | 19              | chr19                        | NC_000019.9      |
| chr19                              | 20               | chr20                         | 20              | chr20                        | NC_000020.10     |
| chr20                              | 21               | chr21                         | 21              | chr21                        | NC_000021.8      |
| chr21                              | 22               | chr22                         | 22              | chr22                        | NC_000022.10     |
| chr22                              | X                | chrX                          | X               | chrX                         | NC_000023.10     |
| chrX                               | Y                | chrY                          | Y               | chrY                         | NC_000024.9      |
| chrY                               | MT               | chrM                          | MT              | chrM                         | NC_012920.1      |
|                                    |                  |                               |                 |                              | NT_113878.1      |
|                                    |                  |                               |                 |                              | NT_113885.1      |
|                                    |                  |                               |                 |                              | NT_113888.1      |
|                                    |                  |                               |                 |                              | NT_113889.1      |
|                                    |                  |                               |                 |                              | NT_113891.2      |
|                                    |                  |                               |                 |                              | NT_113901.1      |
|                                    |                  |                               |                 |                              | NT_113907.1      |
|                                    |                  |                               |                 |                              | NT_113909.1      |
|                                    |                  |                               |                 |                              | NT_113911.1      |
|                                    |                  |                               |                 |                              | NT_113914.1      |
|                                    |                  |                               |                 |                              | NT_113915.1      |
|                                    |                  |                               |                 |                              | NT_113916.2      |
|                                    |                  |                               |                 |                              | NT_113921.2      |
|                                    |                  |                               |                 |                              | NT_113923.1      |
|                                    |                  |                               |                 |                              | NT_113930.1      |
|                                    |                  |                               |                 |                              | NT_113941.1      |
|                                    |                  |                               |                 |                              | NT_113943.1      |
|                                    |                  |                               |                 |                              | NT_113945.1      |
|                                    |                  |                               |                 |                              | NT_113947.1      |
|                                    |                  |                               |                 |                              | NT_113948.1      |
|                                    |                  |                               |                 |                              | NT_113949.1      |
|                                    |                  |                               |                 |                              | NT_113950.2      |
|                                    |                  |                               |                 |                              | NT_113961.1      |
|                                    |                  |                               |                 |                              | NT_167207.1      |
|                                    |                  |                               |                 |                              | NT_167208.1      |
|                                    |                  |                               |                 |                              | NT_167209.1      |
|                                    |                  |                               |                 |                              | NT_167210.1      |
|                                    |                  |                               |                 |                              | NT_167211.1      |
|                                    |                  |                               |                 |                              | NT_167212.1      |
|                                    |                  |                               |                 |                              | NT_167213.1      |
|                                    |                  |                               |                 |                              | NT_167214.1      |
|                                    |                  |                               |                 |                              | NT_167215.1      |
|                                    |                  |                               |                 |                              | NT_167216.1      |
|                                    |                  |                               |                 |                              | NT_167217.1      |
|                                    |                  |                               |                 |                              | NT_167218.1      |
|                                    |                  |                               |                 |                              | NT_167219.1      |
|                                    |                  |                               |                 |                              | NT_167220.1      |
|                                    |                  |                               |                 |                              | NT_167221.1      |
|                                    |                  |                               |                 |                              | NT_167222.1      |
|                                    |                  |                               |                 |                              | NT_167223.1      |
|                                    |                  |                               |                 |                              | NT_167224.1      |
|                                    |                  |                               |                 |                              | NT_167225.1      |
|                                    |                  |                               |                 |                              | NT_167226.1      |
|                                    |                  |                               |                 |                              | NT_167227.1      |
|                                    |                  |                               |                 |                              | NT_167228.1      |
|                                    |                  |                               |                 |                              | NT_167229.1      |
|                                    |                  |                               |                 |                              | NT_167230.1      |
|                                    |                  |                               |                 |                              | NT_167231.1      |
|                                    |                  |                               |                 |                              | NT_167232.1      |
|                                    |                  |                               |                 |                              | NT_167233.1      |
|                                    |                  |                               |                 |                              | NT_167234.1      |
|                                    |                  |                               |                 |                              | NT_167235.1      |
|                                    |                  |                               |                 |                              | NT_167236.1      |
|                                    |                  |                               |                 |                              | NT_167237.1      |
|                                    |                  |                               |                 |                              | NT_167238.1      |
|                                    |                  |                               |                 |                              | NT_167239.1      |
|                                    |                  |                               |                 |                              | NT_167240.1      |
|                                    |                  |                               |                 |                              | NT_167241.1      |
|                                    |                  |                               |                 |                              | NT_167242.1      |
|                                    |                  |                               |                 |                              | NT_167243.1      |
|                                    |                  |                               |                 |                              | NT_167244.1      |
|                                    |                  |                               |                 |                              | NT_167245.1      |
|                                    |                  |                               |                 |                              | NT_167246.1      |
|                                    |                  |                               |                 |                              | NT_167247.1      |
|                                    |                  |                               |                 |                              | NT_167248.1      |
|                                    |                  |                               |                 |                              | NT_167249.1      |
|                                    |                  |                               |                 |                              | NT_167250.1      |
|                                    |                  |                               |                 |                              | NT_167251.1      |
|                                    |                  |                               |                 |                              | NW_003315903.1   |
|                                    |                  |                               |                 |                              | NW_003315904.1   |
|                                    |                  |                               |                 |                              | NW_003315905.1   |
|                                    |                  |                               |                 |                              | NW_003315906.1   |
|                                    |                  |                               |                 |                              | NW_003315907.1   |
|                                    |                  |                               |                 |                              | NW_003315908.1   |
|                                    |                  |                               |                 |                              | NW_003315909.1   |
|                                    |                  |                               |                 |                              | NW_003315910.1   |
|                                    |                  |                               |                 |                              | NW_003315911.1   |
|                                    |                  |                               |                 |                              | NW_003315912.1   |
|                                    |                  |                               |                 |                              | NW_003315913.1   |
|                                    |                  |                               |                 |                              | NW_003315914.1   |
|                                    |                  |                               |                 |                              | NW_003315915.1   |
|                                    |                  |                               |                 |                              | NW_003315916.1   |
|                                    |                  |                               |                 |                              | NW_003315917.2   |
|                                    |                  |                               |                 |                              | NW_003315918.1   |
|                                    |                  |                               |                 |                              | NW_003315919.1   |
|                                    |                  |                               |                 |                              | NW_003315920.1   |
|                                    |                  |                               |                 |                              | NW_003315921.1   |
|                                    |                  |                               |                 |                              | NW_003315922.2   |
|                                    |                  |                               |                 |                              | NW_003315923.1   |
|                                    |                  |                               |                 |                              | NW_003315924.1   |
|                                    |                  |                               |                 |                              | NW_003315925.1   |
|                                    |                  |                               |                 |                              | NW_003315926.1   |
|                                    |                  |                               |                 |                              | NW_003315927.1   |
|                                    |                  |                               |                 |                              | NW_003315928.1   |
|                                    |                  |                               |                 |                              | NW_003315929.1   |
|                                    |                  |                               |                 |                              | NW_003315930.1   |
|                                    |                  |                               |                 |                              | NW_003315931.1   |
|                                    |                  |                               |                 |                              | NW_003315932.1   |
|                                    |                  |                               |                 |                              | NW_003315933.1   |
|                                    |                  |                               |                 |                              | NW_003315934.1   |
|                                    |                  |                               |                 |                              | NW_003315935.1   |
|                                    |                  |                               |                 |                              | NW_003315936.1   |
|                                    |                  |                               |                 |                              | NW_003315937.1   |
|                                    |                  |                               |                 |                              | NW_003315938.1   |
|                                    |                  |                               |                 |                              | NW_003315939.1   |
|                                    |                  |                               |                 |                              | NW_003315940.1   |
|                                    |                  |                               |                 |                              | NW_003315941.1   |
|                                    |                  |                               |                 |                              | NW_003315942.2   |
|                                    |                  |                               |                 |                              | NW_003315943.1   |
|                                    |                  |                               |                 |                              | NW_003315944.1   |
|                                    |                  |                               |                 |                              | NW_003315945.1   |
|                                    |                  |                               |                 |                              | NW_003315946.1   |
|                                    |                  |                               |                 |                              | NW_003315947.1   |
|                                    |                  |                               |                 |                              | NW_003315948.2   |
|                                    |                  |                               |                 |                              | NW_003315949.1   |
|                                    |                  |                               |                 |                              | NW_003315950.2   |
|                                    |                  |                               |                 |                              | NW_003315951.1   |
|                                    |                  |                               |                 |                              | NW_003315952.2   |
|                                    |                  |                               |                 |                              | NW_003315953.1   |
|                                    |                  |                               |                 |                              | NW_003315954.1   |
|                                    |                  |                               |                 |                              | NW_003315955.1   |
|                                    |                  |                               |                 |                              | NW_003315956.1   |
|                                    |                  |                               |                 |                              | NW_003315957.1   |
|                                    |                  |                               |                 |                              | NW_003315958.1   |
|                                    |                  |                               |                 |                              | NW_003315959.1   |
|                                    |                  |                               |                 |                              | NW_003315960.1   |
|                                    |                  |                               |                 |                              | NW_003315961.1   |
|                                    |                  |                               |                 |                              | NW_003315962.1   |
|                                    |                  |                               |                 |                              | NW_003315963.1   |
|                                    |                  |                               |                 |                              | NW_003315964.2   |
|                                    |                  |                               |                 |                              | NW_003315965.1   |
|                                    |                  |                               |                 |                              | NW_003315966.1   |
|                                    |                  |                               |                 |                              | NW_003315967.1   |
|                                    |                  |                               |                 |                              | NW_003315968.1   |
|                                    |                  |                               |                 |                              | NW_003315969.1   |
|                                    |                  |                               |                 |                              | NW_003315970.1   |
|                                    |                  |                               |                 |                              | NW_003315971.2   |
|                                    |                  |                               |                 |                              | NW_003315972.1   |
|                                    |                  |                               |                 |                              | NW_003571030.1   |
|                                    |                  |                               |                 |                              | NW_003571031.1   |
|                                    |                  |                               |                 |                              | NW_003571032.1   |
|                                    |                  |                               |                 |                              | NW_003571033.2   |
|                                    |                  |                               |                 |                              | NW_003571034.1   |
|                                    |                  |                               |                 |                              | NW_003571035.1   |
|                                    |                  |                               |                 |                              | NW_003571036.1   |
|                                    |                  |                               |                 |                              | NW_003571037.1   |
|                                    |                  |                               |                 |                              | NW_003571038.1   |
|                                    |                  |                               |                 |                              | NW_003571039.1   |
|                                    |                  |                               |                 |                              | NW_003571040.1   |
|                                    |                  |                               |                 |                              | NW_003571041.1   |
|                                    |                  |                               |                 |                              | NW_003571042.1   |
|                                    |                  |                               |                 |                              | NW_003571043.1   |
|                                    |                  |                               |                 |                              | NW_003571045.1   |
|                                    |                  |                               |                 |                              | NW_003571046.1   |
|                                    |                  |                               |                 |                              | NW_003571047.1   |
|                                    |                  |                               |                 |                              | NW_003571048.1   |
|                                    |                  |                               |                 |                              | NW_003571049.1   |
|                                    |                  |                               |                 |                              | NW_003571050.1   |
|                                    |                  |                               |                 |                              | NW_003571051.1   |
|                                    |                  |                               |                 |                              | NW_003571052.1   |
|                                    |                  |                               |                 |                              | NW_003571053.2   |
|                                    |                  |                               |                 |                              | NW_003571054.1   |
|                                    |                  |                               |                 |                              | NW_003571055.1   |
|                                    |                  |                               |                 |                              | NW_003571056.1   |
|                                    |                  |                               |                 |                              | NW_003571057.1   |
|                                    |                  |                               |                 |                              | NW_003571058.1   |
|                                    |                  |                               |                 |                              | NW_003571059.1   |
|                                    |                  |                               |                 |                              | NW_003571060.1   |
|                                    |                  |                               |                 |                              | NW_003571061.1   |
|                                    |                  |                               |                 |                              | NW_003571063.2   |
|                                    |                  |                               |                 |                              | NW_003571064.2   |
|                                    |                  |                               |                 |                              | NW_003871055.3   |
|                                    |                  |                               |                 |                              | NW_003871056.3   |
|                                    |                  |                               |                 |                              | NW_003871057.1   |
|                                    |                  |                               |                 |                              | NW_003871058.1   |
|                                    |                  |                               |                 |                              | NW_003871059.1   |
|                                    |                  |                               |                 |                              | NW_003871060.1   |
|                                    |                  |                               |                 |                              | NW_003871061.1   |
|                                    |                  |                               |                 |                              | NW_003871062.1   |
|                                    |                  |                               |                 |                              | NW_003871063.1   |
|                                    |                  |                               |                 |                              | NW_003871064.1   |
|                                    |                  |                               |                 |                              | NW_003871065.1   |
|                                    |                  |                               |                 |                              | NW_003871066.2   |
|                                    |                  |                               |                 |                              | NW_003871067.1   |
|                                    |                  |                               |                 |                              | NW_003871068.1   |
|                                    |                  |                               |                 |                              | NW_003871069.1   |
|                                    |                  |                               |                 |                              | NW_003871070.1   |
|                                    |                  |                               |                 |                              | NW_003871071.1   |
|                                    |                  |                               |                 |                              | NW_003871072.2   |
|                                    |                  |                               |                 |                              | NW_003871073.1   |
|                                    |                  |                               |                 |                              | NW_003871074.1   |
|                                    |                  |                               |                 |                              | NW_003871075.1   |
|                                    |                  |                               |                 |                              | NW_003871076.1   |
|                                    |                  |                               |                 |                              | NW_003871077.1   |
|                                    |                  |                               |                 |                              | NW_003871078.1   |
|                                    |                  |                               |                 |                              | NW_003871079.1   |
|                                    |                  |                               |                 |                              | NW_003871080.1   |
|                                    |                  |                               |                 |                              | NW_003871081.1   |
|                                    |                  |                               |                 |                              | NW_003871082.1   |
|                                    |                  |                               |                 |                              | NW_003871083.2   |
|                                    |                  |                               |                 |                              | NW_003871084.1   |
|                                    |                  |                               |                 |                              | NW_003871085.1   |
|                                    |                  |                               |                 |                              | NW_003871086.1   |
|                                    |                  |                               |                 |                              | NW_003871087.1   |
|                                    |                  |                               |                 |                              | NW_003871088.1   |
|                                    |                  |                               |                 |                              | NW_003871089.1   |
|                                    |                  |                               |                 |                              | NW_003871090.1   |
|                                    |                  |                               |                 |                              | NW_003871091.1   |
|                                    |                  |                               |                 |                              | NW_003871092.1   |
|                                    |                  |                               |                 |                              | NW_003871093.1   |
|                                    |                  |                               |                 |                              | NW_003871094.1   |
|                                    |                  |                               |                 |                              | NW_003871095.1   |
|                                    |                  |                               |                 |                              | NW_003871096.1   |
|                                    |                  |                               |                 |                              | NW_003871098.1   |
|                                    |                  |                               |                 |                              | NW_003871099.1   |
|                                    |                  |                               |                 |                              | NW_003871100.1   |
|                                    |                  |                               |                 |                              | NW_003871101.3   |
|                                    |                  |                               |                 |                              | NW_003871102.1   |
|                                    |                  |                               |                 |                              | NW_003871103.3   |
|                                    |                  |                               |                 |                              | NW_004070863.1   |
|                                    |                  |                               |                 |                              | NW_004070864.2   |
|                                    |                  |                               |                 |                              | NW_004070865.1   |
|                                    |                  |                               |                 |                              | NW_004070866.1   |
|                                    |                  |                               |                 |                              | NW_004070867.1   |
|                                    |                  |                               |                 |                              | NW_004070868.1   |
|                                    |                  |                               |                 |                              | NW_004070869.1   |
|                                    |                  |                               |                 |                              | NW_004070870.1   |
|                                    |                  |                               |                 |                              | NW_004070871.1   |
|                                    |                  |                               |                 |                              | NW_004070872.2   |
|                                    |                  |                               |                 |                              | NW_004070873.1   |
|                                    |                  |                               |                 |                              | NW_004070874.1   |
|                                    |                  |                               |                 |                              | NW_004070875.1   |
|                                    |                  |                               |                 |                              | NW_004070876.1   |
|                                    |                  |                               |                 |                              | NW_004070877.1   |
|                                    |                  |                               |                 |                              | NW_004070878.1   |
|                                    |                  |                               |                 |                              | NW_004070879.1   |
|                                    |                  |                               |                 |                              | NW_004070880.2   |
|                                    |                  |                               |                 |                              | NW_004070881.1   |
|                                    |                  |                               |                 |                              | NW_004070882.1   |
|                                    |                  |                               |                 |                              | NW_004070883.1   |
|                                    |                  |                               |                 |                              | NW_004070884.1   |
|                                    |                  |                               |                 |                              | NW_004070885.1   |
|                                    |                  |                               |                 |                              | NW_004070886.1   |
|                                    |                  |                               |                 |                              | NW_004070887.1   |
|                                    |                  |                               |                 |                              | NW_004070888.1   |
|                                    |                  |                               |                 |                              | NW_004070889.1   |
|                                    |                  |                               |                 |                              | NW_004070890.2   |
|                                    |                  |                               |                 |                              | NW_004070891.1   |
|                                    |                  |                               |                 |                              | NW_004070892.1   |
|                                    |                  |                               |                 |                              | NW_004070893.1   |
|                                    |                  |                               |                 |                              | NW_004166862.1   |
|                                    |                  |                               |                 |                              | NW_004166863.1   |
|                                    |                  |                               |                 |                              | NW_004166864.2   |
|                                    |                  |                               |                 |                              | NW_004166865.1   |
|                                    |                  |                               |                 |                              | NW_004166866.1   |
|                                    |                  |                               |                 |                              | NW_004504299.1   |
|                                    |                  |                               |                 |                              | NW_004504300.1   |
|                                    |                  |                               |                 |                              | NW_004504301.1   |
|                                    |                  |                               |                 |                              | NW_004504302.1   |
|                                    |                  |                               |                 |                              | NW_004504303.2   |
|                                    |                  |                               |                 |                              | NW_004504304.1   |
|                                    |                  |                               |                 |                              | NW_004504305.1   |
|                                    |                  |                               |                 |                              | NW_004775426.1   |
|                                    |                  |                               |                 |                              | NW_004775427.1   |
|                                    |                  |                               |                 |                              | NW_004775428.1   |
|                                    |                  |                               |                 |                              | NW_004775429.1   |
|                                    |                  |                               |                 |                              | NW_004775430.1   |
|                                    |                  |                               |                 |                              | NW_004775431.1   |
|                                    |                  |                               |                 |                              | NW_004775432.1   |
|                                    |                  |                               |                 |                              | NW_004775433.1   |
|                                    |                  |                               |                 |                              | NW_004775434.1   |
|                                    |                  |                               |                 |                              | NW_004775435.1   |

#### Data sources

GATK (build 138)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/dbsnp_138.hg19.vcf.gz
# To produce tabix index file (.tbi)
gunzip dbsnp_138.hg19.vcf.gz
bgzip dbsnp_138.hg19.vcf
tabix dbsnp_138.hg19.vcf.gz
```

NCBI (build 150)

```bash
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz.tbi
```

NCBI (build 150) gatk version

```bash
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/GATK/All_20170710.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/GATK/All_20170710.vcf.gz.tbi
```

NCBI (build151)

```bash
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.tbi
```

NCBI (build151) gatk version

```bash
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/All_20180423.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/All_20180423.vcf.gz.tbi
```

NCBI (build 153)

```bash
wget ftp://ftp.ncbi.nlm.nih.gov:21/snp/latest_release/VCF/GCF_000001405.25.gz
wget ftp://ftp.ncbi.nlm.nih.gov:21/snp/latest_release/VCF/GCF_000001405.25.gz.tbi
```

### GRCh38

| GATK (build 146) | NCBI (build 150) | NCBI (build 150) gatk version | NCBI (build151) | NCBI (build151) gatk version | NCBI (build 153) |
|------------------|------------------|-------------------------------|-----------------|------------------------------|------------------|
| chr1             | 1                | chr1                          | 1               | chr1                         | NC_000001.11     |
| chr2             | 2                | chr2                          | 2               | chr2                         | NC_000002.12     |
| chr3             | 3                | chr3                          | 3               | chr3                         | NC_000003.12     |
| chr4             | 4                | chr4                          | 4               | chr4                         | NC_000004.12     |
| chr5             | 5                | chr5                          | 5               | chr5                         | NC_000005.10     |
| chr6             | 6                | chr6                          | 6               | chr6                         | NC_000006.12     |
| chr7             | 7                | chr7                          | 7               | chr7                         | NC_000007.14     |
| chr8             | 8                | chr8                          | 8               | chr8                         | NC_000008.11     |
| chr9             | 9                | chr9                          | 9               | chr9                         | NC_000009.12     |
| chr10            | 10               | chr10                         | 10              | chr10                        | NC_000010.11     |
| chr11            | 11               | chr11                         | 11              | chr11                        | NC_000011.10     |
| chr12            | 12               | chr12                         | 12              | chr12                        | NC_000012.12     |
| chr13            | 13               | chr13                         | 13              | chr13                        | NC_000013.11     |
| chr14            | 14               | chr14                         | 14              | chr14                        | NC_000014.9      |
| chr15            | 15               | chr15                         | 15              | chr15                        | NC_000015.10     |
| chr16            | 16               | chr16                         | 16              | chr16                        | NC_000016.10     |
| chr17            | 17               | chr17                         | 17              | chr17                        | NC_000017.11     |
| chr18            | 18               | chr18                         | 18              | chr18                        | NC_000018.10     |
| chr19            | 19               | chr19                         | 19              | chr19                        | NC_000019.10     |
| chr20            | 20               | chr20                         | 20              | chr20                        | NC_000020.11     |
| chr21            | 21               | chr21                         | 21              | chr21                        | NC_000021.9      |
| chr22            | 22               | chr22                         | 22              | chr22                        | NC_000022.11     |
| chrX             | X                | chrX                          | X               | chrX                         | NC_000023.11     |
| chrY             | Y                | chrY                          | Y               | chrY                         | NC_000024.10     |
| chrM             | MT               | chrM                          | MT              | chrM                         | NC_012920.1      |
|                  |                  |                               |                 |                              | NT_113793.3      |
|                  |                  |                               |                 |                              | NT_113796.3      |
|                  |                  |                               |                 |                              | NT_113888.1      |
|                  |                  |                               |                 |                              | NT_113889.1      |
|                  |                  |                               |                 |                              | NT_113891.3      |
|                  |                  |                               |                 |                              | NT_113901.1      |
|                  |                  |                               |                 |                              | NT_113930.2      |
|                  |                  |                               |                 |                              | NT_113948.1      |
|                  |                  |                               |                 |                              | NT_113949.2      |
|                  |                  |                               |                 |                              | NT_167208.1      |
|                  |                  |                               |                 |                              | NT_167209.1      |
|                  |                  |                               |                 |                              | NT_167211.2      |
|                  |                  |                               |                 |                              | NT_167213.1      |
|                  |                  |                               |                 |                              | NT_167214.1      |
|                  |                  |                               |                 |                              | NT_167215.1      |
|                  |                  |                               |                 |                              | NT_167218.1      |
|                  |                  |                               |                 |                              | NT_167219.1      |
|                  |                  |                               |                 |                              | NT_167220.1      |
|                  |                  |                               |                 |                              | NT_167244.2      |
|                  |                  |                               |                 |                              | NT_167245.2      |
|                  |                  |                               |                 |                              | NT_167246.2      |
|                  |                  |                               |                 |                              | NT_167247.2      |
|                  |                  |                               |                 |                              | NT_167248.2      |
|                  |                  |                               |                 |                              | NT_167249.2      |
|                  |                  |                               |                 |                              | NT_167250.2      |
|                  |                  |                               |                 |                              | NT_167251.2      |
|                  |                  |                               |                 |                              | NT_187361.1      |
|                  |                  |                               |                 |                              | NT_187362.1      |
|                  |                  |                               |                 |                              | NT_187363.1      |
|                  |                  |                               |                 |                              | NT_187364.1      |
|                  |                  |                               |                 |                              | NT_187365.1      |
|                  |                  |                               |                 |                              | NT_187366.1      |
|                  |                  |                               |                 |                              | NT_187367.1      |
|                  |                  |                               |                 |                              | NT_187368.1      |
|                  |                  |                               |                 |                              | NT_187369.1      |
|                  |                  |                               |                 |                              | NT_187370.1      |
|                  |                  |                               |                 |                              | NT_187371.1      |
|                  |                  |                               |                 |                              | NT_187372.1      |
|                  |                  |                               |                 |                              | NT_187373.1      |
|                  |                  |                               |                 |                              | NT_187374.1      |
|                  |                  |                               |                 |                              | NT_187375.1      |
|                  |                  |                               |                 |                              | NT_187376.1      |
|                  |                  |                               |                 |                              | NT_187377.1      |
|                  |                  |                               |                 |                              | NT_187378.1      |
|                  |                  |                               |                 |                              | NT_187379.1      |
|                  |                  |                               |                 |                              | NT_187380.1      |
|                  |                  |                               |                 |                              | NT_187381.1      |
|                  |                  |                               |                 |                              | NT_187382.1      |
|                  |                  |                               |                 |                              | NT_187383.1      |
|                  |                  |                               |                 |                              | NT_187384.1      |
|                  |                  |                               |                 |                              | NT_187385.1      |
|                  |                  |                               |                 |                              | NT_187386.1      |
|                  |                  |                               |                 |                              | NT_187387.1      |
|                  |                  |                               |                 |                              | NT_187388.1      |
|                  |                  |                               |                 |                              | NT_187389.1      |
|                  |                  |                               |                 |                              | NT_187390.1      |
|                  |                  |                               |                 |                              | NT_187391.1      |
|                  |                  |                               |                 |                              | NT_187392.1      |
|                  |                  |                               |                 |                              | NT_187393.1      |
|                  |                  |                               |                 |                              | NT_187394.1      |
|                  |                  |                               |                 |                              | NT_187395.1      |
|                  |                  |                               |                 |                              | NT_187396.1      |
|                  |                  |                               |                 |                              | NT_187397.1      |
|                  |                  |                               |                 |                              | NT_187398.1      |
|                  |                  |                               |                 |                              | NT_187399.1      |
|                  |                  |                               |                 |                              | NT_187400.1      |
|                  |                  |                               |                 |                              | NT_187401.1      |
|                  |                  |                               |                 |                              | NT_187402.1      |
|                  |                  |                               |                 |                              | NT_187403.1      |
|                  |                  |                               |                 |                              | NT_187405.1      |
|                  |                  |                               |                 |                              | NT_187406.1      |
|                  |                  |                               |                 |                              | NT_187407.1      |
|                  |                  |                               |                 |                              | NT_187408.1      |
|                  |                  |                               |                 |                              | NT_187409.1      |
|                  |                  |                               |                 |                              | NT_187410.1      |
|                  |                  |                               |                 |                              | NT_187411.1      |
|                  |                  |                               |                 |                              | NT_187412.1      |
|                  |                  |                               |                 |                              | NT_187413.1      |
|                  |                  |                               |                 |                              | NT_187414.1      |
|                  |                  |                               |                 |                              | NT_187415.1      |
|                  |                  |                               |                 |                              | NT_187416.1      |
|                  |                  |                               |                 |                              | NT_187417.1      |
|                  |                  |                               |                 |                              | NT_187418.1      |
|                  |                  |                               |                 |                              | NT_187419.1      |
|                  |                  |                               |                 |                              | NT_187420.1      |
|                  |                  |                               |                 |                              | NT_187421.1      |
|                  |                  |                               |                 |                              | NT_187422.1      |
|                  |                  |                               |                 |                              | NT_187423.1      |
|                  |                  |                               |                 |                              | NT_187424.1      |
|                  |                  |                               |                 |                              | NT_187425.1      |
|                  |                  |                               |                 |                              | NT_187426.1      |
|                  |                  |                               |                 |                              | NT_187427.1      |
|                  |                  |                               |                 |                              | NT_187428.1      |
|                  |                  |                               |                 |                              | NT_187429.1      |
|                  |                  |                               |                 |                              | NT_187430.1      |
|                  |                  |                               |                 |                              | NT_187431.1      |
|                  |                  |                               |                 |                              | NT_187432.1      |
|                  |                  |                               |                 |                              | NT_187433.1      |
|                  |                  |                               |                 |                              | NT_187434.1      |
|                  |                  |                               |                 |                              | NT_187435.1      |
|                  |                  |                               |                 |                              | NT_187436.1      |
|                  |                  |                               |                 |                              | NT_187437.1      |
|                  |                  |                               |                 |                              | NT_187438.1      |
|                  |                  |                               |                 |                              | NT_187439.1      |
|                  |                  |                               |                 |                              | NT_187440.1      |
|                  |                  |                               |                 |                              | NT_187441.1      |
|                  |                  |                               |                 |                              | NT_187443.1      |
|                  |                  |                               |                 |                              | NT_187444.1      |
|                  |                  |                               |                 |                              | NT_187446.1      |
|                  |                  |                               |                 |                              | NT_187447.1      |
|                  |                  |                               |                 |                              | NT_187448.1      |
|                  |                  |                               |                 |                              | NT_187449.1      |
|                  |                  |                               |                 |                              | NT_187450.1      |
|                  |                  |                               |                 |                              | NT_187451.1      |
|                  |                  |                               |                 |                              | NT_187452.1      |
|                  |                  |                               |                 |                              | NT_187453.1      |
|                  |                  |                               |                 |                              | NT_187454.1      |
|                  |                  |                               |                 |                              | NT_187455.1      |
|                  |                  |                               |                 |                              | NT_187456.1      |
|                  |                  |                               |                 |                              | NT_187457.1      |
|                  |                  |                               |                 |                              | NT_187458.1      |
|                  |                  |                               |                 |                              | NT_187461.1      |
|                  |                  |                               |                 |                              | NT_187464.1      |
|                  |                  |                               |                 |                              | NT_187465.1      |
|                  |                  |                               |                 |                              | NT_187466.1      |
|                  |                  |                               |                 |                              | NT_187467.1      |
|                  |                  |                               |                 |                              | NT_187469.1      |
|                  |                  |                               |                 |                              | NT_187470.1      |
|                  |                  |                               |                 |                              | NT_187471.1      |
|                  |                  |                               |                 |                              | NT_187473.1      |
|                  |                  |                               |                 |                              | NT_187476.1      |
|                  |                  |                               |                 |                              | NT_187478.1      |
|                  |                  |                               |                 |                              | NT_187480.1      |
|                  |                  |                               |                 |                              | NT_187481.1      |
|                  |                  |                               |                 |                              | NT_187482.1      |
|                  |                  |                               |                 |                              | NT_187484.1      |
|                  |                  |                               |                 |                              | NT_187485.1      |
|                  |                  |                               |                 |                              | NT_187486.1      |
|                  |                  |                               |                 |                              | NT_187488.1      |
|                  |                  |                               |                 |                              | NT_187490.1      |
|                  |                  |                               |                 |                              | NT_187493.1      |
|                  |                  |                               |                 |                              | NT_187495.1      |
|                  |                  |                               |                 |                              | NT_187496.1      |
|                  |                  |                               |                 |                              | NT_187497.1      |
|                  |                  |                               |                 |                              | NT_187498.1      |
|                  |                  |                               |                 |                              | NT_187499.1      |
|                  |                  |                               |                 |                              | NT_187500.1      |
|                  |                  |                               |                 |                              | NT_187501.1      |
|                  |                  |                               |                 |                              | NT_187502.1      |
|                  |                  |                               |                 |                              | NT_187503.1      |
|                  |                  |                               |                 |                              | NT_187504.1      |
|                  |                  |                               |                 |                              | NT_187505.1      |
|                  |                  |                               |                 |                              | NT_187506.1      |
|                  |                  |                               |                 |                              | NT_187508.1      |
|                  |                  |                               |                 |                              | NT_187509.1      |
|                  |                  |                               |                 |                              | NT_187510.1      |
|                  |                  |                               |                 |                              | NT_187511.1      |
|                  |                  |                               |                 |                              | NT_187512.1      |
|                  |                  |                               |                 |                              | NT_187513.1      |
|                  |                  |                               |                 |                              | NT_187514.1      |
|                  |                  |                               |                 |                              | NT_187515.1      |
|                  |                  |                               |                 |                              | NT_187516.1      |
|                  |                  |                               |                 |                              | NT_187517.1      |
|                  |                  |                               |                 |                              | NT_187518.1      |
|                  |                  |                               |                 |                              | NT_187519.1      |
|                  |                  |                               |                 |                              | NT_187520.1      |
|                  |                  |                               |                 |                              | NT_187521.1      |
|                  |                  |                               |                 |                              | NT_187522.1      |
|                  |                  |                               |                 |                              | NT_187523.1      |
|                  |                  |                               |                 |                              | NT_187524.1      |
|                  |                  |                               |                 |                              | NT_187525.1      |
|                  |                  |                               |                 |                              | NT_187526.1      |
|                  |                  |                               |                 |                              | NT_187527.1      |
|                  |                  |                               |                 |                              | NT_187528.1      |
|                  |                  |                               |                 |                              | NT_187529.1      |
|                  |                  |                               |                 |                              | NT_187530.1      |
|                  |                  |                               |                 |                              | NT_187531.1      |
|                  |                  |                               |                 |                              | NT_187532.1      |
|                  |                  |                               |                 |                              | NT_187533.1      |
|                  |                  |                               |                 |                              | NT_187534.1      |
|                  |                  |                               |                 |                              | NT_187535.1      |
|                  |                  |                               |                 |                              | NT_187536.1      |
|                  |                  |                               |                 |                              | NT_187537.1      |
|                  |                  |                               |                 |                              | NT_187538.1      |
|                  |                  |                               |                 |                              | NT_187539.1      |
|                  |                  |                               |                 |                              | NT_187540.1      |
|                  |                  |                               |                 |                              | NT_187541.1      |
|                  |                  |                               |                 |                              | NT_187542.1      |
|                  |                  |                               |                 |                              | NT_187543.1      |
|                  |                  |                               |                 |                              | NT_187544.1      |
|                  |                  |                               |                 |                              | NT_187545.1      |
|                  |                  |                               |                 |                              | NT_187546.1      |
|                  |                  |                               |                 |                              | NT_187547.1      |
|                  |                  |                               |                 |                              | NT_187548.1      |
|                  |                  |                               |                 |                              | NT_187549.1      |
|                  |                  |                               |                 |                              | NT_187550.1      |
|                  |                  |                               |                 |                              | NT_187551.1      |
|                  |                  |                               |                 |                              | NT_187552.1      |
|                  |                  |                               |                 |                              | NT_187553.1      |
|                  |                  |                               |                 |                              | NT_187554.1      |
|                  |                  |                               |                 |                              | NT_187555.1      |
|                  |                  |                               |                 |                              | NT_187556.1      |
|                  |                  |                               |                 |                              | NT_187557.1      |
|                  |                  |                               |                 |                              | NT_187558.1      |
|                  |                  |                               |                 |                              | NT_187559.1      |
|                  |                  |                               |                 |                              | NT_187560.1      |
|                  |                  |                               |                 |                              | NT_187561.1      |
|                  |                  |                               |                 |                              | NT_187562.1      |
|                  |                  |                               |                 |                              | NT_187563.1      |
|                  |                  |                               |                 |                              | NT_187564.1      |
|                  |                  |                               |                 |                              | NT_187565.1      |
|                  |                  |                               |                 |                              | NT_187566.1      |
|                  |                  |                               |                 |                              | NT_187567.1      |
|                  |                  |                               |                 |                              | NT_187568.1      |
|                  |                  |                               |                 |                              | NT_187569.1      |
|                  |                  |                               |                 |                              | NT_187570.1      |
|                  |                  |                               |                 |                              | NT_187571.1      |
|                  |                  |                               |                 |                              | NT_187572.1      |
|                  |                  |                               |                 |                              | NT_187573.1      |
|                  |                  |                               |                 |                              | NT_187574.1      |
|                  |                  |                               |                 |                              | NT_187575.1      |
|                  |                  |                               |                 |                              | NT_187576.1      |
|                  |                  |                               |                 |                              | NT_187577.1      |
|                  |                  |                               |                 |                              | NT_187578.1      |
|                  |                  |                               |                 |                              | NT_187579.1      |
|                  |                  |                               |                 |                              | NT_187580.1      |
|                  |                  |                               |                 |                              | NT_187581.1      |
|                  |                  |                               |                 |                              | NT_187582.1      |
|                  |                  |                               |                 |                              | NT_187583.1      |
|                  |                  |                               |                 |                              | NT_187584.1      |
|                  |                  |                               |                 |                              | NT_187585.1      |
|                  |                  |                               |                 |                              | NT_187586.1      |
|                  |                  |                               |                 |                              | NT_187587.1      |
|                  |                  |                               |                 |                              | NT_187588.1      |
|                  |                  |                               |                 |                              | NT_187589.1      |
|                  |                  |                               |                 |                              | NT_187590.1      |
|                  |                  |                               |                 |                              | NT_187591.1      |
|                  |                  |                               |                 |                              | NT_187592.1      |
|                  |                  |                               |                 |                              | NT_187593.1      |
|                  |                  |                               |                 |                              | NT_187594.1      |
|                  |                  |                               |                 |                              | NT_187595.1      |
|                  |                  |                               |                 |                              | NT_187596.1      |
|                  |                  |                               |                 |                              | NT_187597.1      |
|                  |                  |                               |                 |                              | NT_187598.1      |
|                  |                  |                               |                 |                              | NT_187599.1      |
|                  |                  |                               |                 |                              | NT_187600.1      |
|                  |                  |                               |                 |                              | NT_187601.1      |
|                  |                  |                               |                 |                              | NT_187602.1      |
|                  |                  |                               |                 |                              | NT_187603.1      |
|                  |                  |                               |                 |                              | NT_187604.1      |
|                  |                  |                               |                 |                              | NT_187605.1      |
|                  |                  |                               |                 |                              | NT_187606.1      |
|                  |                  |                               |                 |                              | NT_187607.1      |
|                  |                  |                               |                 |                              | NT_187608.1      |
|                  |                  |                               |                 |                              | NT_187609.1      |
|                  |                  |                               |                 |                              | NT_187610.1      |
|                  |                  |                               |                 |                              | NT_187611.1      |
|                  |                  |                               |                 |                              | NT_187612.1      |
|                  |                  |                               |                 |                              | NT_187613.1      |
|                  |                  |                               |                 |                              | NT_187614.1      |
|                  |                  |                               |                 |                              | NT_187615.1      |
|                  |                  |                               |                 |                              | NT_187616.1      |
|                  |                  |                               |                 |                              | NT_187617.1      |
|                  |                  |                               |                 |                              | NT_187618.1      |
|                  |                  |                               |                 |                              | NT_187619.1      |
|                  |                  |                               |                 |                              | NT_187620.1      |
|                  |                  |                               |                 |                              | NT_187621.1      |
|                  |                  |                               |                 |                              | NT_187622.1      |
|                  |                  |                               |                 |                              | NT_187623.1      |
|                  |                  |                               |                 |                              | NT_187624.1      |
|                  |                  |                               |                 |                              | NT_187625.1      |
|                  |                  |                               |                 |                              | NT_187626.1      |
|                  |                  |                               |                 |                              | NT_187627.1      |
|                  |                  |                               |                 |                              | NT_187628.1      |
|                  |                  |                               |                 |                              | NT_187629.1      |
|                  |                  |                               |                 |                              | NT_187630.1      |
|                  |                  |                               |                 |                              | NT_187631.1      |
|                  |                  |                               |                 |                              | NT_187632.1      |
|                  |                  |                               |                 |                              | NT_187633.1      |
|                  |                  |                               |                 |                              | NT_187634.1      |
|                  |                  |                               |                 |                              | NT_187635.1      |
|                  |                  |                               |                 |                              | NT_187636.1      |
|                  |                  |                               |                 |                              | NT_187637.1      |
|                  |                  |                               |                 |                              | NT_187638.1      |
|                  |                  |                               |                 |                              | NT_187639.1      |
|                  |                  |                               |                 |                              | NT_187640.1      |
|                  |                  |                               |                 |                              | NT_187641.1      |
|                  |                  |                               |                 |                              | NT_187642.1      |
|                  |                  |                               |                 |                              | NT_187643.1      |
|                  |                  |                               |                 |                              | NT_187644.1      |
|                  |                  |                               |                 |                              | NT_187645.1      |
|                  |                  |                               |                 |                              | NT_187646.1      |
|                  |                  |                               |                 |                              | NT_187647.1      |
|                  |                  |                               |                 |                              | NT_187648.1      |
|                  |                  |                               |                 |                              | NT_187649.1      |
|                  |                  |                               |                 |                              | NT_187650.1      |
|                  |                  |                               |                 |                              | NT_187651.1      |
|                  |                  |                               |                 |                              | NT_187652.1      |
|                  |                  |                               |                 |                              | NT_187653.1      |
|                  |                  |                               |                 |                              | NT_187654.1      |
|                  |                  |                               |                 |                              | NT_187655.1      |
|                  |                  |                               |                 |                              | NT_187656.1      |
|                  |                  |                               |                 |                              | NT_187657.1      |
|                  |                  |                               |                 |                              | NT_187658.1      |
|                  |                  |                               |                 |                              | NT_187659.1      |
|                  |                  |                               |                 |                              | NT_187660.1      |
|                  |                  |                               |                 |                              | NT_187661.1      |
|                  |                  |                               |                 |                              | NT_187662.1      |
|                  |                  |                               |                 |                              | NT_187663.1      |
|                  |                  |                               |                 |                              | NT_187664.1      |
|                  |                  |                               |                 |                              | NT_187665.1      |
|                  |                  |                               |                 |                              | NT_187666.1      |
|                  |                  |                               |                 |                              | NT_187667.1      |
|                  |                  |                               |                 |                              | NT_187668.1      |
|                  |                  |                               |                 |                              | NT_187669.1      |
|                  |                  |                               |                 |                              | NT_187670.1      |
|                  |                  |                               |                 |                              | NT_187671.1      |
|                  |                  |                               |                 |                              | NT_187672.1      |
|                  |                  |                               |                 |                              | NT_187673.1      |
|                  |                  |                               |                 |                              | NT_187674.1      |
|                  |                  |                               |                 |                              | NT_187675.1      |
|                  |                  |                               |                 |                              | NT_187676.1      |
|                  |                  |                               |                 |                              | NT_187677.1      |
|                  |                  |                               |                 |                              | NT_187678.1      |
|                  |                  |                               |                 |                              | NT_187679.1      |
|                  |                  |                               |                 |                              | NT_187680.1      |
|                  |                  |                               |                 |                              | NT_187681.1      |
|                  |                  |                               |                 |                              | NT_187682.1      |
|                  |                  |                               |                 |                              | NT_187683.1      |
|                  |                  |                               |                 |                              | NT_187684.1      |
|                  |                  |                               |                 |                              | NT_187685.1      |
|                  |                  |                               |                 |                              | NT_187686.1      |
|                  |                  |                               |                 |                              | NT_187687.1      |
|                  |                  |                               |                 |                              | NT_187688.1      |
|                  |                  |                               |                 |                              | NT_187689.1      |
|                  |                  |                               |                 |                              | NT_187690.1      |
|                  |                  |                               |                 |                              | NT_187691.1      |
|                  |                  |                               |                 |                              | NT_187692.1      |
|                  |                  |                               |                 |                              | NT_187693.1      |
|                  |                  |                               |                 |                              | NW_003315905.1   |
|                  |                  |                               |                 |                              | NW_003315906.1   |
|                  |                  |                               |                 |                              | NW_003315907.2   |
|                  |                  |                               |                 |                              | NW_003315908.1   |
|                  |                  |                               |                 |                              | NW_003315909.1   |
|                  |                  |                               |                 |                              | NW_003315913.1   |
|                  |                  |                               |                 |                              | NW_003315914.1   |
|                  |                  |                               |                 |                              | NW_003315915.1   |
|                  |                  |                               |                 |                              | NW_003315917.2   |
|                  |                  |                               |                 |                              | NW_003315918.1   |
|                  |                  |                               |                 |                              | NW_003315919.1   |
|                  |                  |                               |                 |                              | NW_003315920.1   |
|                  |                  |                               |                 |                              | NW_003315921.1   |
|                  |                  |                               |                 |                              | NW_003315922.2   |
|                  |                  |                               |                 |                              | NW_003315928.1   |
|                  |                  |                               |                 |                              | NW_003315929.1   |
|                  |                  |                               |                 |                              | NW_003315930.1   |
|                  |                  |                               |                 |                              | NW_003315931.1   |
|                  |                  |                               |                 |                              | NW_003315934.1   |
|                  |                  |                               |                 |                              | NW_003315935.1   |
|                  |                  |                               |                 |                              | NW_003315936.1   |
|                  |                  |                               |                 |                              | NW_003315938.1   |
|                  |                  |                               |                 |                              | NW_003315939.2   |
|                  |                  |                               |                 |                              | NW_003315940.1   |
|                  |                  |                               |                 |                              | NW_003315941.1   |
|                  |                  |                               |                 |                              | NW_003315942.2   |
|                  |                  |                               |                 |                              | NW_003315943.1   |
|                  |                  |                               |                 |                              | NW_003315944.2   |
|                  |                  |                               |                 |                              | NW_003315945.1   |
|                  |                  |                               |                 |                              | NW_003315946.1   |
|                  |                  |                               |                 |                              | NW_003315952.3   |
|                  |                  |                               |                 |                              | NW_003315953.2   |
|                  |                  |                               |                 |                              | NW_003315954.1   |
|                  |                  |                               |                 |                              | NW_003315955.1   |
|                  |                  |                               |                 |                              | NW_003315956.1   |
|                  |                  |                               |                 |                              | NW_003315957.1   |
|                  |                  |                               |                 |                              | NW_003315958.1   |
|                  |                  |                               |                 |                              | NW_003315959.1   |
|                  |                  |                               |                 |                              | NW_003315960.1   |
|                  |                  |                               |                 |                              | NW_003315961.1   |
|                  |                  |                               |                 |                              | NW_003315962.1   |
|                  |                  |                               |                 |                              | NW_003315963.1   |
|                  |                  |                               |                 |                              | NW_003315964.2   |
|                  |                  |                               |                 |                              | NW_003315965.1   |
|                  |                  |                               |                 |                              | NW_003315966.2   |
|                  |                  |                               |                 |                              | NW_003315967.2   |
|                  |                  |                               |                 |                              | NW_003315968.2   |
|                  |                  |                               |                 |                              | NW_003315969.2   |
|                  |                  |                               |                 |                              | NW_003315970.2   |
|                  |                  |                               |                 |                              | NW_003315971.2   |
|                  |                  |                               |                 |                              | NW_003315972.2   |
|                  |                  |                               |                 |                              | NW_003571033.2   |
|                  |                  |                               |                 |                              | NW_003571036.1   |
|                  |                  |                               |                 |                              | NW_003571049.1   |
|                  |                  |                               |                 |                              | NW_003571050.1   |
|                  |                  |                               |                 |                              | NW_003571054.1   |
|                  |                  |                               |                 |                              | NW_003571055.2   |
|                  |                  |                               |                 |                              | NW_003571056.2   |
|                  |                  |                               |                 |                              | NW_003571057.2   |
|                  |                  |                               |                 |                              | NW_003571058.2   |
|                  |                  |                               |                 |                              | NW_003571059.2   |
|                  |                  |                               |                 |                              | NW_003571060.1   |
|                  |                  |                               |                 |                              | NW_003571061.2   |
|                  |                  |                               |                 |                              | NW_003871060.2   |
|                  |                  |                               |                 |                              | NW_003871073.1   |
|                  |                  |                               |                 |                              | NW_003871074.1   |
|                  |                  |                               |                 |                              | NW_003871091.1   |
|                  |                  |                               |                 |                              | NW_003871092.1   |
|                  |                  |                               |                 |                              | NW_003871093.1   |
|                  |                  |                               |                 |                              | NW_004166862.2   |
|                  |                  |                               |                 |                              | NW_004504305.1   |
|                  |                  |                               |                 |                              | NW_009646194.1   |
|                  |                  |                               |                 |                              | NW_009646195.1   |
|                  |                  |                               |                 |                              | NW_009646196.1   |
|                  |                  |                               |                 |                              | NW_009646197.1   |
|                  |                  |                               |                 |                              | NW_009646198.1   |
|                  |                  |                               |                 |                              | NW_009646199.1   |
|                  |                  |                               |                 |                              | NW_009646200.1   |
|                  |                  |                               |                 |                              | NW_009646201.1   |
|                  |                  |                               |                 |                              | NW_009646202.1   |
|                  |                  |                               |                 |                              | NW_009646203.1   |
|                  |                  |                               |                 |                              | NW_009646204.1   |
|                  |                  |                               |                 |                              | NW_009646205.1   |
|                  |                  |                               |                 |                              | NW_009646206.1   |
|                  |                  |                               |                 |                              | NW_009646207.1   |
|                  |                  |                               |                 |                              | NW_009646208.1   |
|                  |                  |                               |                 |                              | NW_009646209.1   |
|                  |                  |                               |                 |                              | NW_011332687.1   |
|                  |                  |                               |                 |                              | NW_011332688.1   |
|                  |                  |                               |                 |                              | NW_011332689.1   |
|                  |                  |                               |                 |                              | NW_011332690.1   |
|                  |                  |                               |                 |                              | NW_011332691.1   |
|                  |                  |                               |                 |                              | NW_011332692.1   |
|                  |                  |                               |                 |                              | NW_011332693.1   |
|                  |                  |                               |                 |                              | NW_011332694.1   |
|                  |                  |                               |                 |                              | NW_011332695.1   |
|                  |                  |                               |                 |                              | NW_011332696.1   |
|                  |                  |                               |                 |                              | NW_011332697.1   |
|                  |                  |                               |                 |                              | NW_011332698.1   |
|                  |                  |                               |                 |                              | NW_011332699.1   |
|                  |                  |                               |                 |                              | NW_011332700.1   |
|                  |                  |                               |                 |                              | NW_011332701.1   |
|                  |                  |                               |                 |                              | NW_012132914.1   |
|                  |                  |                               |                 |                              | NW_012132915.1   |
|                  |                  |                               |                 |                              | NW_012132916.1   |
|                  |                  |                               |                 |                              | NW_012132917.1   |
|                  |                  |                               |                 |                              | NW_012132918.1   |
|                  |                  |                               |                 |                              | NW_012132919.1   |
|                  |                  |                               |                 |                              | NW_012132920.1   |
|                  |                  |                               |                 |                              | NW_012132921.1   |
|                  |                  |                               |                 |                              | NW_013171799.1   |
|                  |                  |                               |                 |                              | NW_013171800.1   |
|                  |                  |                               |                 |                              | NW_013171801.1   |
|                  |                  |                               |                 |                              | NW_013171802.1   |
|                  |                  |                               |                 |                              | NW_013171803.1   |
|                  |                  |                               |                 |                              | NW_013171804.1   |
|                  |                  |                               |                 |                              | NW_013171805.1   |
|                  |                  |                               |                 |                              | NW_013171806.1   |
|                  |                  |                               |                 |                              | NW_013171807.1   |
|                  |                  |                               |                 |                              | NW_013171808.1   |
|                  |                  |                               |                 |                              | NW_013171809.1   |
|                  |                  |                               |                 |                              | NW_013171810.1   |
|                  |                  |                               |                 |                              | NW_013171811.1   |
|                  |                  |                               |                 |                              | NW_013171812.1   |
|                  |                  |                               |                 |                              | NW_013171813.1   |
|                  |                  |                               |                 |                              | NW_013171814.1   |
|                  |                  |                               |                 |                              | NW_014040925.1   |
|                  |                  |                               |                 |                              | NW_014040926.1   |
|                  |                  |                               |                 |                              | NW_014040927.1   |
|                  |                  |                               |                 |                              | NW_014040928.1   |
|                  |                  |                               |                 |                              | NW_014040929.1   |
|                  |                  |                               |                 |                              | NW_014040930.1   |
|                  |                  |                               |                 |                              | NW_014040931.1   |
|                  |                  |                               |                 |                              | NW_015148966.1   |
|                  |                  |                               |                 |                              | NW_015148967.1   |
|                  |                  |                               |                 |                              | NW_015148968.1   |
|                  |                  |                               |                 |                              | NW_015148969.1   |
|                  |                  |                               |                 |                              | NW_015495298.1   |
|                  |                  |                               |                 |                              | NW_015495299.1   |
|                  |                  |                               |                 |                              | NW_015495300.1   |
|                  |                  |                               |                 |                              | NW_015495301.1   |
|                  |                  |                               |                 |                              | NW_016107297.1   |
|                  |                  |                               |                 |                              | NW_016107298.1   |
|                  |                  |                               |                 |                              | NW_016107299.1   |
|                  |                  |                               |                 |                              | NW_016107300.1   |
|                  |                  |                               |                 |                              | NW_016107301.1   |
|                  |                  |                               |                 |                              | NW_016107302.1   |
|                  |                  |                               |                 |                              | NW_016107303.1   |
|                  |                  |                               |                 |                              | NW_016107304.1   |
|                  |                  |                               |                 |                              | NW_016107305.1   |
|                  |                  |                               |                 |                              | NW_016107306.1   |
|                  |                  |                               |                 |                              | NW_016107307.1   |
|                  |                  |                               |                 |                              | NW_016107308.1   |
|                  |                  |                               |                 |                              | NW_016107309.1   |
|                  |                  |                               |                 |                              | NW_016107310.1   |
|                  |                  |                               |                 |                              | NW_016107311.1   |
|                  |                  |                               |                 |                              | NW_016107312.1   |
|                  |                  |                               |                 |                              | NW_016107313.1   |
|                  |                  |                               |                 |                              | NW_016107314.1   |
|                  |                  |                               |                 |                              | NW_017363813.1   |
|                  |                  |                               |                 |                              | NW_017363814.1   |
|                  |                  |                               |                 |                              | NW_017363815.1   |
|                  |                  |                               |                 |                              | NW_017363816.1   |
|                  |                  |                               |                 |                              | NW_017363817.1   |
|                  |                  |                               |                 |                              | NW_017363818.1   |
|                  |                  |                               |                 |                              | NW_017363819.1   |
|                  |                  |                               |                 |                              | NW_017363820.1   |
|                  |                  |                               |                 |                              | NW_017852928.1   |
|                  |                  |                               |                 |                              | NW_017852929.1   |
|                  |                  |                               |                 |                              | NW_017852930.1   |
|                  |                  |                               |                 |                              | NW_017852931.1   |
|                  |                  |                               |                 |                              | NW_017852932.1   |
|                  |                  |                               |                 |                              | NW_017852933.1   |
|                  |                  |                               |                 |                              | NW_018654706.1   |
|                  |                  |                               |                 |                              | NW_018654707.1   |
|                  |                  |                               |                 |                              | NW_018654708.1   |
|                  |                  |                               |                 |                              | NW_018654709.1   |
|                  |                  |                               |                 |                              | NW_018654710.1   |
|                  |                  |                               |                 |                              | NW_018654711.1   |
|                  |                  |                               |                 |                              | NW_018654712.1   |
|                  |                  |                               |                 |                              | NW_018654713.1   |
|                  |                  |                               |                 |                              | NW_018654714.1   |
|                  |                  |                               |                 |                              | NW_018654715.1   |
|                  |                  |                               |                 |                              | NW_018654716.1   |
|                  |                  |                               |                 |                              | NW_018654717.1   |
|                  |                  |                               |                 |                              | NW_018654718.1   |
|                  |                  |                               |                 |                              | NW_018654719.1   |
|                  |                  |                               |                 |                              | NW_018654720.1   |
|                  |                  |                               |                 |                              | NW_018654721.1   |
|                  |                  |                               |                 |                              | NW_018654722.1   |
|                  |                  |                               |                 |                              | NW_018654723.1   |
|                  |                  |                               |                 |                              | NW_018654724.1   |
|                  |                  |                               |                 |                              | NW_018654725.1   |
|                  |                  |                               |                 |                              | NW_018654726.1   |
|                  |                  |                               |                 |                              | NW_019805487.1   |
|                  |                  |                               |                 |                              | NW_019805488.1   |
|                  |                  |                               |                 |                              | NW_019805489.1   |
|                  |                  |                               |                 |                              | NW_019805490.1   |
|                  |                  |                               |                 |                              | NW_019805491.1   |
|                  |                  |                               |                 |                              | NW_019805492.1   |
|                  |                  |                               |                 |                              | NW_019805493.1   |
|                  |                  |                               |                 |                              | NW_019805494.1   |
|                  |                  |                               |                 |                              | NW_019805495.1   |
|                  |                  |                               |                 |                              | NW_019805496.1   |
|                  |                  |                               |                 |                              | NW_019805497.1   |
|                  |                  |                               |                 |                              | NW_019805498.1   |
|                  |                  |                               |                 |                              | NW_019805499.1   |
|                  |                  |                               |                 |                              | NW_019805500.1   |
|                  |                  |                               |                 |                              | NW_019805501.1   |
|                  |                  |                               |                 |                              | NW_019805502.1   |
|                  |                  |                               |                 |                              | NW_019805503.1   |

#### Data sources

GATK (build 146)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/dbsnp_146.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi
```

NCBI (build 150)

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/All_20170710.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/All_20170710.vcf.gz.tbi
```

NCBI (build 150) gatk version

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/All_20170710.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/All_20170710.vcf.gz.tbi
```

NCBI (build 151)

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz.tbi
```

NCBI (build 151) gatk version

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/All_20180418.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/All_20180418.vcf.gz.tbi
```

NCBI (build 153)

```bash
wget ftp://ftp.ncbi.nlm.nih.gov:21/snp/latest_release/VCF/GCF_000001405.38.gz
wget ftp://ftp.ncbi.nlm.nih.gov:21/snp/latest_release/VCF/GCF_000001405.38.gz.tbi
```
