####################
### The goals of this script are:
### To filter the SNP data using VCFtools.
### To extract one SNP per locus.
### To change data format from vcf to geno (012) for sNMF analyses.
### Checked for functionality on April 20th 2021.

### Good resource for VCFtools: https://www.ddocent.com/filtering/

## Selecting our data (name based on ipyrad parameters):
for run in atlas_gr_n96_ms04 colletti_gr_n27_ms04 essingtonii_gr_n40_ms04 inornatus_gr_n218_ms04 leonhardii_gr_n133_ms04 pantherinus_gr_n73_ms04 schomburgkii_gr_n101_ms04 taeniolatus_gr_n11_ms04
do

## Maximum missing data per SNP position:
ms=0.4

## Creating directories:
mkdir ~/Dropbox/Science/MYPAPERS_ongoing/spheno_ddRAD/VCFtools/ms${ms}
mkdir ~/Dropbox/Science/MYPAPERS_ongoing/spheno_ddRAD/VCFtools/ms${ms}/${run}
cd ~/Dropbox/Science/MYPAPERS_ongoing/spheno_ddRAD/VCFtools/ms${ms}/${run}

## Filtering SNPs by MAF and bialelic (diploid) condition.
vcftools --vcf ~/Dropbox/Science/MYPAPERS_ongoing/spheno_ddRAD/ipyrad_outfiles/ms${ms}/${run}_outfiles/${run}.vcf --recode --out snps --maf 0.05 --max-alleles 2

## Renaming resulting vcf file:
mv snps.recode.vcf snps.vcf

## Now extracting one SNP per locus:
cat snps.vcf | grep "#" > usnps.vcf ## Extract the headers (whose line starts with #), save in different file.
cat snps.vcf | grep -v "#" | sort -u -k1,1 >> usnps.vcf ## Sort by RAD locus name and extract unique SNP per locus; concatenate with that file created.

## Saving in 012 format:
vcftools --vcf usnps.vcf --out usnps --012

## Replacing -1 with 9 to indicate missing sites:
cat usnps.012 | sed -e 's/-1/9/g' usnps.012 > usnps_tmp.012 ## Replace and save in temporary file.

## Lastly, renaming temporary file as a final 012 file:
mv usnps_tmp.012 usnps.012

## Closing loop:
done

## Renaming folders:
cd ..
rename 's/_n.+//' *

## End of script.
