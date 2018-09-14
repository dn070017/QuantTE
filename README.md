 # QuantTE
 ## About
 QuantTE is an analsis pipeline that utilize the <b>RNA-Seq</b> data and <b>RepeatMasker table</b> to quantify the repeat element. There are two main stages for QuantTE: <b>extraction</b> and <b>quantification</b>. The first stage used the RepeatMasker table and genome sequences to extract the transposable element in interest while the second stage used Kallisto to quantify the abundance of transposable element.
 ## Basic Command
 - Extract stage:
 ```shell
 python3 QuantTE/main.py --extract --input input.json
 ```
 - Quant stage:
 ```shell
 python3 QuantTE/main.py --kallisto --input input.json
 ```
 - Run Complete Analysis
 ```shell
 python3 QuantTE/main.py --extract --kallisto --input input.json
 ```
 ## Preparation of Input.json 
 - Complete analysis (Single-End Reads)
 ```json
 {
    "output_dir": "path of output directory",
    "genome_fasta": "path of genomic sequence",
    "transcript_fasta": "path of transcriptomic sequence",
    "repeat_masker_dir": "directory of repeat masker tables",
    "read_label": ["condition.1", "condition.2"],
    "read_fastq": ["condition.1.fastq", "condition.2.fastq"],
    "TE_fasta": "(output) path of transposon element sequence"
 }
 ```
 - Complete analysis (Paired-End Reads)
 ```json
 {
    "output_dir": "path of output directory",
    "genome_fasta": "path of genomic sequence",
    "transcript_fasta": "path of transcriptomic sequence",
    "repeat_masker_dir": "directory of repeat masker tables",
    "read_label": ["condition.1", "condition.2"],
    "read_fastq": ["condition.1.R1.fastq condition.1.R2.fastq", "condition.2.R1.fastq condition.2.R2.fastq"],
    "TE_fasta": "(output) path of transposon element sequence"
 }
 ```


