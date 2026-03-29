## Overview

This script calculates linkage disequilibrium (LD) between a target SNP and all variants within a specified genomic window using PLINK.

### Software download

Download: <https://www.cog-genomics.org/plink/>

### Installation

``` bash
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip

unzip plink_linux_x86_64_20250819.zip
```

### Add PLINK to PATH

This allows you to run PLINK from any directory without specifying the full path

1.  Open you `~/.bashrc`

2.  Export thr plink path to the end of the file

    `export PATH="/path/to/plink/download:$PATH"`

3.  Reload your `.bashrc`

    `source ~/.bashrc`

4.  Verify PLINK is accessible

    `plink –help`

### Usage

```         
bash extractLD_bfile.sh <plink_path> <bfile_prefix> <chr> <pos_start> <pos_end> <top_snp> <window_kb> <out_file>
```

### Example (5kb window)

For a SNP at position 955672 on chromosome Pf3D7_09_v3:

``` bash
bash extractLD_bfile.sh /home/bofosu/Software/plink ../CLEAN_RESULTS/BFILE/bed_file Pf3D7_09_v3 950672 960672 Pf3D7_09_v3:955672 5 ld_output
```

![](images/clipboard-653814727.png)
