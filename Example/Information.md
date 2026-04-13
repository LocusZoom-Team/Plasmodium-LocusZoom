### Information Regarding the Test Run for LocusZoom Plasmodium Falciparum

### LD Script

Chromosome: 7

SNP: 296804

Code:

```{bash}
./plink --bfile /home/user/Plasmodium-LocusZoom/LD_Bed/bed_file --allow-extra-chr --chr Pf3D7_07_v3 --from-bp 250000 --to-bp 350000 --r2 --ld-window 100000000 --ld-window-kb 500 --ld-window-r2 0 --ld-snp Pf3D7_07_v3:296804 --keep-allele-order --out /home/user/plinker/testrunLD
```
### LocusZoom Run

```{bash}
Rscript test.R example.assoc.log testrunLD.ld Pf3D7_07_v3:296804 Example.genes Pf3D7_07_v3:309298 10000 SNP_296804_Association_With_Mono /home/user/Plasmodium-LocusZoom/Example/output.jpeg
```