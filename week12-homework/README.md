# Week 12 Homework

## Step 1.

### Step 1B. 

I used the provided code, saved as ```parse_kraken.py```

```
for KRAKEN in 83 86 88 89 90 93 94 97; do ./parse_kraken.py ./metagenomics_data/step0_givendata/KRAKEN/SRR4921${KRAKEN}.kraken SRR4921${KRAKEN}; done
```

### Step 1C. 

Rank ktImportText on the samples

```
ktImportText -q *_krona.txt
```

Throughout the first week, _Enterococcus faecalis_ is the most prevalent group of bacteria, with the other major groups being _Bacillales_ (mostly including different _staphyloccocus_ strains) and _Actinobacteria_. On the first day the amount of _Bacillales_ and _Actinobacteria_ are approximately the same (~18% of total) with the rest being _Enterococcus faecalis_. Then the _Actinobacteria_ group collapses entirely and _Bacillales_ decreases, with _Enterococcus faecalis_ taking up 92% of the samples. Over the rest of the week, _Bacillales_ gradually increases again, and in the last 2 days, _Actinobacteria_ reappears and then rapidly increases to 27% on the last day.