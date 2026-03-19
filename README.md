# COI Recovery Pipeline using MitoGeneExtractor and ORF-based Filtering

## Overview

This pipeline provides a reproducible workflow to recover mitochondrial COI sequences from paired-end sequencing data (e.g., UCE, transcriptomes, or shotgun data), followed by biologically informed filtering based on open reading frame (ORF) integrity.

The workflow includes:

1. Preparation of input FASTQ files
2. Retrieval and curation of COI reference sequences (GenBank)
3. COI recovery using MitoGeneExtractor
4. Post-processing of consensus sequences
5. ORF validation using mitochondrial genetic code
6. Selection of the best sequences per sample
7. Standardized renaming of FASTA files

---

## Directory Structure

mitogene_COI/
├── fastq/        # Paired-end reads (*READ1.fastq / *READ2.fastq)
├── refs/         # Reference COI sequences
├── results/      # Output from MitoGeneExtractor
└── README.md

---

## Step 0 — Prepare Input FASTQ Files

Ensure all reads are in plain FASTQ format (not compressed).

If your files are compressed (.fastq.gz), decompress them:

```
gunzip *.fastq.gz
```

Expected naming pattern:

sample-READ1.fastq
sample-READ2.fastq

---

## Step 1 — Build COI Reference Dataset (GenBank)

Retrieve COI sequences for your target taxon (e.g., Theraphosidae) from GenBank.

Recommended procedure:

* Query example: Theraphosidae[Organism] AND COI
* Download sequences in FASTA format
* Remove duplicates and very short sequences
* Retain representative diversity across taxa

Save as:

refs/COI_reference.fasta

---

## Step 2 — COI Recovery with MitoGeneExtractor

```
for r1 in fastq/*READ1.fastq
do
    sample=$(basename "$r1" -READ1.fastq)

    MitoGeneExtractor-v1.9.6 \
        -q fastq/${sample}-READ1.fastq \
        -q fastq/${sample}-READ2.fastq \
        -p refs/COI_reference.fasta \
        -o results/${sample}_align \
        -c results/${sample}_consensus \
        -C 5 \
        -t 0.2 \
        -r 1 \
        -n 0

done
```

---

## Step 3 — Collect Consensus Sequences

```
mkdir step1_consensus
cp results/*_consensus*.fas step1_consensus/
cd step1_consensus
```

---

## Step 4 — Remove Empty FASTA Files

```
for f in *.fas; do
    seq=$(grep -v ">" "$f" | tr -d '\n')
    if [ -z "$seq" ]; then
        echo "$f"
    fi
done > empty_files.txt

mkdir empty_backup
mv $(cat empty_files.txt) empty_backup/
```

---

## Step 5 — ORF Detection

Create script: find_best_orf.sh

```
#!/bin/bash

mkdir -p orf_results

for f in *.fas; do
    base=$(basename $f .fas)

    best_stop=999
    best_len=0

    for strand in fwd rev; do

        if [ "$strand" = "rev" ]; then
            revseq $f tmp.fas >/dev/null 2>&1
            seqfile="tmp.fas"
        else
            seqfile=$f
        fi

        for frame in 1 2 3; do

            transeq -table 5 -frame $frame -sequence $seqfile -outseq tmp.pep >/dev/null 2>&1

            stops=$(grep -o "*" tmp.pep | wc -l)
            len=$(grep -v ">" tmp.pep | tr -d '\n' | wc -c)

            if [ $stops -lt $best_stop ] || ([ $stops -eq $best_stop ] && [ $len -gt $best_len ]); then
                best_stop=$stops
                best_len=$len
            fi

        done
    done

    echo -e "$base\tstops=$best_stop\tlen=$best_len" >> orf_results/summary.tsv

done
```

Run:

```
chmod +x find_best_orf.sh
./find_best_orf.sh
```

---

## Step 6 — Select Best Sequences per Sample

```
awk '{
split($2,a,"=");
split($3,b,"=");
print $1"\t"a[2]"\t"b[2]
}' orf_results/summary.tsv > clean.tsv

awk '{
split($1,a,"_");
id=a[1];
print id"\t"$0
}' clean.tsv > clean_with_id.tsv

sort -k1,1 -k3,3n -k4,4nr clean_with_id.tsv > sorted.tsv

awk '{
id=$1;
stops=$3;
len=$4;

if (stops <= 1 && len >= 200) {
    count[id]++;
    if (count[id] <= 3) {
        print $0
    }
}
}' sorted.tsv > final_selection.tsv
```

---

## Step 7 — Extract Selected FASTA

```
mkdir FINAL_COI

for name in $(cut -f2 final_selection.tsv); do
    cp ${name}.fas FINAL_COI/
done

cd FINAL_COI
```

---

## Step 8 — Rename FASTA Files

```
for f in *_consensus*.fas; do

    base=$(basename "$f" .fas)

    tmp=$(echo "$base" | sed 's/_consensus/_/')
    tmp=$(echo "$tmp" | sed 's/\.1_1//')

    acc=$(echo "$tmp" | awk -F "_" '{print $NF}')
    prefix=$(echo "$tmp" | sed "s/_$acc$//")

    newname="${prefix}_${acc}.fas"

    sed "1s/.*/>${prefix}_${acc}/" "$f" > "$newname"

done
```

---

## Step 9 — Alignment

```
cat *.fas > all_COI.fas
mafft all_COI.fas > aligned.fas
```

---

## Software and Attribution

COI sequences were initially recovered using MitoGeneExtractor (v1.9.6).

This pipeline adds ORF-based validation and selection to improve sequence quality.

Users should cite the original MitoGeneExtractor publication when using this workflow.

---

## Author

Tiago Belintani **Brave the Sun**


