# flagger.smk
## Nov. 3, 2022
This is a Snakemake copy of [**WDL Flagger**](https://github.com/mobinasri/flagger).
## Step 1. Prereqs
```
ln -s ./runlocal .
ln -s ./runsnake .
cp -r config/config.yaml .
```
The directory should look like this:
```
runsnake
runlocal
config
├── config.yaml
└── manifest.tab
```
## Step 2. Inputs
Use the copied config/config.yaml as config with no need to modify.
### config/manifest.tab
```bash
sample	hap1	hap2	fofn	tech
HG00733	path/to/hap1.fasta	path/to/hap2.fasta	path/to/sample_fastq.fofn	hifi
HG00733	path/to/hap1.fasta	path/to/hap2.fasta	path/to/sample_fastq.fofn	ont
```
## Step 3. Start the analysis!
```bash
./runsnake 30
```