# Running alignment script

Usually I can "pseudo parallelize" it as follows:

```
align_task() {
    /broad/thechenlab/Dawn/MyScripts/align_single_fastq_to_reference_George_TRACE_pythonPileup.sh \
        $1 \
        04_RUNX1.fasta
}
N=30
cp /broad/thechenlab/Dawn/RNAreference/04_RUNX1.fasta .
bowtie2-build 04_RUNX1.fasta 04_RUNX1
(
    for FILE in DC575*_R1_001.fastq.gz; do
        ((i = i % N))
        ((i++ == 0)) && wait
        align_task "$FILE" &
    done
)
```


