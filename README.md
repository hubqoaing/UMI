Before running, run ```HTSeq``` first:

```bash
samtools view -o b2_bar.sam                          \
    accepted_hits.bam &&                             \
python htseq-count                                   \
    -s no -f sam -a 10 -o b2_bar.umi.sort.gene.sam   \
    b2_bar.sam                                       \
    refGene.gtf                                      \
    >b2_bar.deseq.xls
```

Then run this script:

```python
python 04_tpm_hubq.py  -g list_gene.xls -e list_ercc.xls b2_bar.umi.sort.gene.sam
```