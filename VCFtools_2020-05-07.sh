####
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### May 2020.

### The goals of this script are:
### To filter the SNP data by MAP using VCFTools.
### To change data formatting to get them ready for STRUCTURE-type of analyses (e.g., sNMF, ADMIXTURE).

### Good resource for VCFtools: https://www.ddocent.com/filtering/

## Creating a folder to save our outputs into:
mkdir /home/ivan/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/VCFtools

## Selecting our data (name based on ipyrad parameters):
for s in inornatus_gr_c90_ni281_mi070 #inornatus_gr_c90_ni281_mi028 
do

## Creating a folder and changing directory to that new folder:
mkdir ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/VCFtools/${s}
cd ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/VCFtools/${s}

## Setting MAC:
for mac in 2 3 4
do

## Creating a folder and changing directory to that new folder:
mkdir ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/VCFtools/${s}/mac${mac}
cd ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/VCFtools/${s}/mac${mac}

## Let's first exclude the outgroups:
vcftools --vcf ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/ipyrad_outfiles/${s}_outfiles/${s}.vcf --remove ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/sample_information/outgroups.indv --recode --out ${s}.noouts

## Let's check levels of individual missing data for ingroup samples:
vcftools --vcf ${s}.noouts.recode.vcf --missing-indv

## Individual missing data (mid) allowed:
for mid in 0.8 0.85
do 

## Creating a folder and changing directory to that new folder:
mkdir ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/VCFtools/${s}/mac${mac}/mid${mid}
cd ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/VCFtools/${s}/mac${mac}/mid${mid}

## Let's create a list of individuals with more than mid (miss_ind) missing sites. We'll use this info to exclude them later.
## We'll use a little trick: An external script to circumvent the fact that mawk won't take ${}.
## The script contains (a command corresponding to) the following: mawk '$5 > ${mid}' out.imiss | cut -f1 > lowDP.indv
bash ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/scripts/mawk_indv_miss_data_mid${mid}.sh

## Add to this list of samples to be excluded four additional problematic samples:
cat ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/sample_information/exclude.txt >> lowDP.indv

## SNP missing data (ms) allowed:
for ms in 0.3 0.33 0.35 0.4
do

## Creating a folder and changing directory to that new folder:
mkdir ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/VCFtools/${s}/mac${mac}/mid${mid}/ms${ms}
cd ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/VCFtools/${s}/mac${mac}/mid${mid}/ms${ms}

## Filtering by MAF/MAC condition and minimum presence data per SNP:
vcftools --vcf ../../*.noouts.recode.vcf --max-missing ${ms} --recode --out ms${ms}.filtered1 --mac ${mac}

## Filter individuals with more than our threshold of missing data:
vcftools --vcf ms${ms}.filtered1.recode.vcf --remove ../lowDP.indv --recode --out ms${ms}.filtered2

## Now extracting one random SNP per locus:
cat ms${ms}.filtered2.recode.vcf | grep "#" > ms${ms}.usnps.vcf # extract the headers (whose line starts with #), save in different file
cat ms${ms}.filtered2.recode.vcf | grep -v "#" | sort -u -k1,1 >> ms${ms}.usnps.vcf # sort by 'chromosome' name and extract unique, concatenate to that file created.

## Saving in 012 format:
vcftools --vcf ms${ms}.usnps.vcf --out ms${ms}.usnps --012

## Replacing -1 with 9 to indicate missing sites:
cat ms${ms}.usnps.012 | sed -e 's/-1/9/g' ms${ms}.usnps.012 > ms${ms}.usnps_tmp.012 # replace and save in temporary file

## Lastly, renaming temporary file as a final 012 file:
mv ms${ms}.usnps_tmp.012 ms${ms}.usnps.012 

## We can also save in plink format to input into ADMIXTURE:
#cat ms${ms}.usnps.vcf | grep ^# > ms${ms}.usnps.plink.vcf 
#cat ms${ms}.usnps.vcf | grep -v ^# | sed 's/^/a/' >> ms${ms}.usnps.plink.vcf
#~/plink2 --vcf ms${ms}.usnps.plink.vcf --make-bed --out ms${ms}.usnps --allow-extra-chr --double-id

## ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0:
#awk '{$1=0;print $0}' ms${ms}.usnps.bim > ms${ms}.usnps.bim.tmp
#mv ms${ms}.usnps.bim.tmp ms${ms}.usnps.bim

## Merging log files:
cat ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/VCFtools/${s}/mac${mac}/*.noouts.log *.filtered1.log *.filtered2.log *.usnps.log > VCFtools_ms${ms}.screen

## Removing temporary files and some other files that we won't need:
rm *.filtered1.recode.vcf
rm *.filtered2.recode.vcf
rm *.pos
rm *.filtered1.log
rm *.filtered2.log
rm *.usnps.log

## Closing loop (ms values):
done

## Closing loop (mid values):
done

## Removing more files we won't need:
rm ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/VCFtools/${s}/mac${mac}/*.log
rm ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/VCFtools/${s}/mac${mac}/*.noouts.recode.vcf

## Closing loop (mac values):
done

## Closing loop (s values):
done

## End of script.
