#!/bin/bash

### set scripts directory ###
export scripts_path='/home/dongke/TF_RNA/other_species/scripts'

### usage document ###
usage () {
	echo "Paremeter illustration are in the following:"
	echo "		-a	TFs Fasta file path		The vcf file name. You can claim the vcf file or gzvcf file folowing this paremeter."
	echo "		-d	DBD region file		you should define the population1 name follow this paremeter, which should lead to the sample list file and will be used for output file name."
	echo "		-i	IDRs score file		you should define the population2 name follow this paremeter, which should lead to the sample list file and will be used for output file name."
	echo "		-f	TFs family index file. "
	echo "		-s	TFs family statistic file. "
	1>&2; exit 1;
}

### get paremeter ###
while getopts "a:d:i:f:s:" opt
do
	case $opt in
	a)
		fasta_file=$OPTARG;;
	d)
		dbd_file=$OPTARG;;
	i)
		idr_file=$OPTARG;;
	f)
		family_index=$OPTARG;;
	s)
		family_sta=$OPTARG;;
	?)
		echo "there is wrong paremter ${OPTARG}" 
		usage ; exit 1 ;;
	esac
done


### verify paremeter ###
if [ -z $fasta_file ] || [ -z $dbd_file ] || [ -z $idr_file ] || [ -z $family_index ] || [ -z $family_sta ] ; then usage ; exit 1 
else
	echo "TFs Fasta file : $fasta_file"
	echo "DBD region file : $dbd_file"
	echo "IDRs score file : $idr_file"
	echo "TFs family index file : $family_index"
	echo "TFs family statistic file: $family_sta"
fi


## R/K rich basic patches scan ##
perl ${scripts_path}/bin/01.R.K.scan.pl $fasta_file

perl ${scripts_path}/bin/02.basic.patches.filtered.by.idr.pl RK.raw.basic.patches.info $idr_file


## DBD ##
# interproscan
#bash /home/dongke/BioSoftware/InterProScan/interproscan-5.59-91.0/interproscan.sh -i $fasta_file

# extract DBD by EntryID
perl ${scripts_path}/bin/03.extract.DBD.region.info.pl $dbd_file

# adjust DBD region by IDRs
perl ${scripts_path}/bin/04.DBD.region.adjust.by.IDRs.pl ${dbd_file}.DBD $idr_file


## Statistic of family frequency and distance ##
perl ${scripts_path}/bin/05.sta.Family.frequency.and.distance.pl RK.raw.basic.patches.info.filtered.by.IDRs ${dbd_file}.DBD.adjusted.by.IDRs.new.region $family_sta


## figure ##
# merge dbd idr correlation
python ${scripts_path}/figure/01.merge.idr.dbd.ARM_Cor.plot.py -f $family_index -i $idr_file -d ${dbd_file}.DBD.adjusted.by.IDRs.new.region -a $fasta_file

# generate family counts and ratio barplot
python ${scripts_path}/figure/02.TF.Family.counts.ratio.barplot.py RK.raw.basic.patches.info.filtered.by.IDRs.distance.to.dbd.lt.20.family.sta