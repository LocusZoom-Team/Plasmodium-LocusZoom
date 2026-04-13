[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5154379.svg)](https://doi.org/10.5281/zenodo.5154379)

# Authors
Original LocusZoom Authors: Tanya J Major and Riku Takei

_Plasmodium falciparum_ Authors: Kashif Bokhari, Belinda Ofosu, and Gabriel Wegman

## Originally LocusZoom-like plots for human GWAS, Now for _Plasmodoium Falciparum_ GWAS
This package was built from the original LocusZoom plotting tool: https://github.com/Geeketics/LocusZooms

Original Paper On LocusZoom: https://doi.org/10.1093/bioinformatics/btq419

The purpose of this tool is to allow users to generate LocusZoom plots from _Plasmodium falciparum_ GWAS data. This package does not generate plots for non-Plasmodium data.

## Make LocusZoom-like plots for _Plasmodium flaciparum_ using MalariaGen and our LD Matrix Generation File

********
The GWAS data used for generating LD matrices can be found at https://www.malariagen.net/resource/34/

The Phenotypic data needed for logistic association metrics can be found in the Phenotypes folder.

********

This package allows the user to create regional Manhattan plots from p-values, log(p-values), or log(Bayes Factors) with points coloured according to LD and genes annotated beneath. The LD input can be generated from the users own data (e.g. for a non-reference population). The package comes with a number of reference files for gene annotation, but is not limited to the use of these files.

This script creates an R function to create LocusZoom-like plots. Three example input files are included for test purposes, along with an example .jpg output.

This script has one package dependency: `scales`

### LocusZoom Requirements

  - Logistic Association: A file of PLINK association results (only the "CHR", "SNP", "BP", and "P" columns are essential)
  - Linkage Disequillibrium: A file of the LD between the SNP to be labelled (top-hit / SNP of interest) and the SNPs included in the PLINK results file
    - this file MUST have a column called "SNP_B" (containing a list of all the SNPs in the results file) and a column called "R2" (containing the R^2 LD value of each SNP). The SNP names MUST match the names in the SNP column of the results file.
  - Gene List: A file of the genes within the region for use in the annotation step. This file must have five columns, "Gene", "Chrom", "Start", "End", and "Coding". The `Pfalciparum_gene_list.txt` file in the `genelist` folder can be used for this file.

### Generating Requirements

  - Logistic Assocation: Files in the `Phenotypes` folder can be used for your LocusZoom runs. If you wish to generate your own, see requirements above.
  - Linkage Disequillibrium: You will generate your own LD file for each run of LocusZoom based on what SNP you are measuring the LD of. Use the file `extractLD_bfile.sh` within the `functions` to help generate the LD file.
  - Gene List: A Plasmodium falciparum gene list `Pfalciparum_gene_list.txt` is provided in the `genelist` folder. If you wish to generate your own from a P. falciparum (.gff) file, use `GetGeneList.py`, also present in the `genelist` folder.

### Generating Your Own LD via Our LD Script

You must have PLINK installed to generate your own LD file. PLINK can be installed via:

```{bash}
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip
unzip plink_linux_x86_64_20230116.zip
```
You must have three bed files generated from your genotype VCF results. This can be done through PLINK. Select the desired chromosome from `Pfalciparum_gene_list.txt`. This must be an exact match to the name present in the row/column (i.e. Pf3D7_10_v3 selects chromosome 10). Select your desired base pair range via the --from-bp and --to-bp arguments. Set a range of 50-1000 for the --ld-window-kb argument. To select your lead SNP (--ld-snp), make sure that it is coming from the same chromosome, and that it exists within the basepair window specified.

Example terminal code: 

```{bash}
./plink --bfile /home/user/Plasmodium-LocusZoom/LD_Bed/bed_file --allow-extra-chr --chr Pf3D7_07_v3 --from-bp 250000 --to-bp 350000 --r2 --ld-window 100000000 --ld-window-kb 500 --ld-window-r2 0 --ld-snp Pf3D7_07_v3:296804 --keep-allele-order --out /home/user/LD_files/testrunLD
```
### For example LocusZoom-like plot generation, see `Information.md` in the `Example` folder