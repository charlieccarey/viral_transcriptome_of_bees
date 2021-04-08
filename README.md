# viral_transcriptome_of_bees
*de novo* Transcriptome assembly of bee viruses.

Bees herein refer to either honey bee *Apis mellifera* or *Andrena spp.*

This repository is in support of the manuscript:

"Metatranscriptome Analysis of Sympatric Bee Species Identifies Bee Virus Variants and a New Virus, Andrena-Associated Bee Virus-1." _Viruses_. 2021 Feb; 13(2): 291.

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7917660/
doi: 10.3390/v13020291
PMCID: PMC7917660
PMID: 33673324

13 samples (described in the accompanying publication):
- were sequenced (Illumina HiSeq 4000 and MiSeq).
- The sequence was processed to count and remove reads associated with certain species. 
- The remaining reads (enriched for virus reads) were assembled with Trinity (trinityrnaseq). 
- Contigs were given putative identities by performing blast and diamond searches against various databases.
