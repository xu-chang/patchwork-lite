# patchwork-lite

Modified from the R project "Patchwork" on http://patchwork.r-forge.r-project.org/. 
patchwork-lite can read VCF results from newer versions of Samtools and BcfTools (Tested on v1.3 on both) without the need for mpileup files to lower the disk-space requirements.

* Compared to the original tool, I made change to patchwork.plot.r and patchwork.alleledata.r, added a python script mpile2alleles.py for reading VCF files.
* Rd documentation files are not changed.

## Usage：
Prepare the VCF file like this:
```bash
samtools mpileup -t 'SP,INFO/AD,INFO/ADF,INFO/ADR' -uvf hg19.fa -r chr2:100000-110000 test.bam | \
  bcftools call --ploidy GRCh37 -mv | \
  bcftools filter -g 5 -G 5 -e 'TYPE!="snp" || AD[0]<5 || DP4[1]<2 || DP4[2]<2' -O z -o ./test.vcf.gz
```
so that we have only SNV lines and DP4 field is in the file.

And then run plot function like this:
```R
patchwork.plot(Tumor_bam, Tumor_vcf, Normal_bam, Normal_vcf, Alpha=alpha, SD=sd)
```

