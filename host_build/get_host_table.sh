#!/bin/bash
# based on https://www.biostars.org/p/306380/
shopt -s globstar


while [[ $# -gt 0 ]]
   do

	TAXID="$1"
	host_table="assembly_summary_refseq.txt"
	patric_table="patric_host_summary.txt"
	if [[ -z "$TAXID" ]]; then
		echo "Usage: build_host_genomes.sh NcbiTaxId"
		exit
	else
		if [[ ! -e $host_table ]]; then
			echo "downloading table"
			wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
		fi
		echo "creating patric table"
		grep '^#' $host_table > $patric_table
		awk -v TAXID="$TAXID" \
			'BEGIN{FS="\t"}{if($6==TAXID){print $0; exit}}' \
			$host_table >> $patric_table
		reference=$(awk -v TAXID="$TAXID" \
			'BEGIN{FS="\t"}{if($6==TAXID){print $20; exit}}' \
			$host_table)
		echo $reference
		base_url=$( printf '%s' "$reference" | awk 'BEGIN{OFS=FS="/"}{print $0,$NF"_genomic"}' )
		mkdir -p $TAXID
		echo "hello"
		base_file=$( echo $base_url | sed 's|.*/||')
		base_file="${TAXID}/${base_file}"
		fna_gz="${base_file}_genomic.fna.gz"
		fna_file="${base_file}_genomic.fna"
		gff_file="${base_file}_genomic.gff"
		gff_gz="${base_file}_genomic.gff.gz"
		if [[ ! -e $fna_gz ]] && [[ ! -e $fna_file ]]; then
			echo $reference | awk 'BEGIN{OFS=FS="/"}{print $0,$NF"_genomic.fna.gz"}' | xargs -n1 wget -O $fna_gz
			gunzip -c $fna_gz > $fna_file
		fi
		if [[ ! -e $gff_gz ]] && [[ ! -e $gff_file ]]; then
			echo $reference | awk 'BEGIN{OFS=FS="/"}{print $0,$NF"_genomic.gff.gz"}' | xargs -n1 wget $O $gff_gz
			gunzip -c $gff_gz > $gff_file
		fi
		hisat_files=( ./$base_file*.ht2 )
		if [[ ${#hisat_files[@]} == 0 ]]; then
			echo "building hisat2 files for ${TAXID}"
			hisat2-build $fna_file $base_file
		fi
	fi

shift
done
