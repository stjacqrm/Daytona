#!/bin/bash 

#----------
# Defaults
#----------

quantitative='FALSE'
from_tar=''
from_dir=''
output_dir='daytona_output'
daytona_variation='clearlabs'
pangolin_docker='docker://staphb/pangolin:4.3.1-pdata-1.25.1'
pangolin_version='v4.3.1'
pdata='v1.25.1'
sotc='S:L452R,S:E484K'
vadr_docker='docker://staphb/vadr:1.3'
samtools_docker='docker://staphb/samtools:1.12'
nextclade_docker='docker://nextstrain/nextclade:0.14.1'

#-------
# Usage
#-------

display_usage() {
  echo -e "\nUsage: run_daytona.sh -t from_tar -r from_dir -o output_dir -v 'daytona_variation' -d 'pangolin_docker' -V 'pangolin_version' -D 'pdata' -s 'sotc' -a 'vadr_docker' -S 'samtools_docker' -n 'nextclade_docker' \n"
  echo -e "  -t  if running from tar file."
  echo -e "      Default: none "
  echo -e "  -r  if running from previously analyzed run."
  echo -e "      Default: none "
  echo -e "  -o  name of Daytona output directory."
  echo -e "      Default: 'daytona_output'"
  echo -e "  -v  name of Daytona Pipeline Variant."
  echo -e "      Default: 'clearlabs'"
  echo -e "  -d  Pangolin docker container."
  echo -e "      Default: 'docker://staphb/pangolin:4.3.1-pdata-1.25.1'"
  echo -e "  -V  Pangolin Version."
  echo -e "      Default: 'v4.3.1'"
  echo -e "  -D  Pangolearn Version."
  echo -e "      Default: 'v1.25.1'"
  echo -e "  -s  SOTCs to evaluate for."
  echo -e "      Default: 'S:L452R,S:E484K'"
  echo -e "  -a  VADR docker container."
  echo -e "      Default: 'docker://staphb/vadr:1.3'"
  echo -e "  -S  Samtools docker container."
  echo -e "      Default: 'docker://staphb/samtools:1.12'"
  echo -e "  -n  Nextclade docker container."
  echo -e "      Default: 'docker://nextstrain/nextclade:0.14.1'"
  echo -e "  -h  Displays this help message\n"
  echo -e "run_daytona.sh"
  echo -e "This wrapper is used to tidy and analyse Clear Labs data "
  echo -e "using a custom cleaner and the Daytona pipeline.\n"
  exit
}


#---------
# Options
#---------

while getopts "t:r:o:v:d:V:D:s:a:S:n:h" opt; do
  case $opt in
    t) from_tar=${OPTARG%/};;
    r) from_dir=${OPTARG%/};;
    o) output_dir=${OPTARG%/};;
    v) daytona_variation=${OPTARG%/};;
    d) pangolin_docker=${OPTARG%/};;
    V) pangolin_version=${OPTARG%/};;
    s) sotc=${OPTARG%/};;
    D) pdata=${OPTARG%/};;
    a) vadr_docker=${OPTARG%/};;
    S) samtools_docker=${OPTARG%/};;
    n) nextclade_docker=${OPTARG%/};;
    h) display_usage;;
   \?) #unrecognized option - show help
      echo -e \\n"Option -$OPTARG not allowed."
      display_usage;;
  esac
done

if [ "$from_tar" = "tar" ]; then
  #-------------------
  # Running File Tidy
  #-------------------

  # making directory name from tarfile name
  filename=$(ls *.tar)
  dir_name=${filename%.all*}
  echo "$dir_name"

  echo 'Extracting Files'
  tar xf *.tar

  echo 'Fixing Filenames'
  for fname in *; do
    name="${fname%\.*}"
    extension="${fname#$name}"
    newname="${name//./_}"
    newfname="$newname""$extension"
    if [ "$fname" != "$newfname" ]; then
      #echo mv "$fname" "$newfname"
      mv "$fname" "$newfname"
    fi
  done

  for fname in *\ *; do 
    mv "$fname" "${fname// /_}"; 
  done

  echo 'Sorting files into appropriate directories'
  mkdir fasta && mv *.fasta fasta
  mkdir fastq && mv *.fastq fastq
  mkdir bam && mv *.bam bam
  mkdir index_bam && mv *.bai index_bam
  mkdir vcf && mv *.vcf vcf
  mkdir metrics && mv *.csv metrics

  echo 'Moving directories into data'
  mkdir data
  mv bam data
  mv fasta data
  mv fastq data
  mv index_bam data
  mv vcf data

  echo 'Removing tar file'
  rm *.tar

  mkdir $dir_name
  mv * $dir_name
  cd $dir_name/metrics && awk 'FNR==1 && NR!=1{next;}{print}' *.csv > multi_metrics.csv
  awk 'BEGIN {FS=","; OFS="\t"} {$1=$1; print}' multi_metrics.csv > multi_metrics1.tsv
  sed -e 's/ /_/g' multi_metrics1.tsv > multi_metrics.tsv
  rm multi_metrics.csv && rm multi_metrics1.tsv
  mv multi_metrics.tsv .. && cd .. && mv multi_metrics.tsv metrics && cd ..

  multi_mets=$dir_name/metrics/multi_metrics.tsv

  echo 'Done Tidying'

  #---------------------------------------
  # Generating params_clearlabs.yaml file
  #---------------------------------------

  echo 'Generating params_clearlabs.yaml file'

  get_abs_filename() {
    # $1 : relative filename
    echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
  }

  fastq="/data/fastq"
  alignment="/data/bam"
  fasta="/data/fasta"


  base=$(get_abs_filename $dir_name)
  outbase=$(get_abs_filename $output_dir)

  input="\"${base}${fastq}\""
  bams="\"${base}${alignment}\""
  assemblies="\"${base}${fasta}\""
  output="\"${outbase}\""
  echo "$output"
  sotc="\"$sotc\""

  pangolin_docker="\"$pangolin_docker\""
  pangolin_version="\"$pangolin_version\""
  pangolin_data_version="\"$pdata\""

  vadr_docker="\"$vadr_docker\""
  samtools_docker="\"$samtools_docker\""
  nc_d=$nextclade_docker
  nextclade_docker="\"$nextclade_docker\""

  touch params_clearlabs.yaml
  echo "# The parameters "input", "bams", "assemblies" and "output" are the absolute paths of the folders. Do not include the "/" at the end of the paths." >> params_clearlabs.yaml
  echo "input : "$input | tee -a params_clearlabs.yaml
  echo "bams : "$bams | tee -a params_clearlabs.yaml
  echo "assemblies : "$assemblies | tee -a params_clearlabs.yaml
  echo "output : "$output | tee -a params_clearlabs.yaml
  echo ""| tee -a params_clearlabs.yaml
  echo "# comma separated list of SOTCs to screen, default as "S:L452R,S:E484K""| tee -a params_clearlabs.yaml
  echo "sotc : "$sotc | tee -a params_clearlabs.yaml
  echo ""| tee -a params_clearlabs.yaml
  echo "#Do not change the parameter settings below:" | tee -a params_clearlabs.yaml
  echo "pangolin_docker : "$pangolin_docker | tee -a params_clearlabs.yaml
  echo "pangolin_version : "$pangolin_version | tee -a params_clearlabs.yaml
  echo "pangolin_data_version : "$pangolin_data_version | tee -a params_clearlabs.yaml

  echo ""| tee -a params_clearlabs.yaml
  echo "vadr_docker : "$vadr_docker | tee -a params_clearlabs.yaml
  echo "samtools_docker : "$samtools_docker | tee -a params_clearlabs.yaml
  echo "nextclade_docker : "$nextclade_docker | tee -a params_clearlabs.yaml
 
  #---------------------------
  # Running Daytona Pipeline
  #---------------------------

  if [ "$daytona_variation" = "clearlabs" ]; then
    ##### run clearlabs pipeline
    echo "run clearlabs pipeline"
    nextflow run ~/Applications/daytona/flaq_sc2_clearlabs2.nf -params-file params_clearlabs.yaml

    sort $output_dir/*/report.txt | uniq > $output_dir/sum_report.txt
    sed -i '/sampleID\treference/d' $output_dir/sum_report.txt
    sed -i '1i sampleID\treference\tstart\tend\tnum_clean_reads\tnum_mapped_reads\tpercent_mapped_clean_reads\tcov_bases_mapped\tpercent_genome_cov_map\tmean_depth\tmean_base_qual\tmean_map_qual\tassembly_length\tnumN\tpercent_ref_genome_cov\tVADR_flag\tQC_flag\tpangolin_version\tlineage\tSOTC' $output_dir/sum_report.txt

    awk 'BEGIN {FS=OFS="\t"} {print $1,$8,$10,$14,$15,$17,$18,$19}' $output_dir/sum_report.txt > $output_dir/mid_file
    paste -d'\t' $multi_mets $output_dir/mid_file | awk 'BEGIN {FS=OFS="\t"} {print $3,$4,$5,$6,$7,$8,$9,$10,$2}' > $output_dir/daytona_report.tsv
    rm $output_dir/mid_file

    cat $output_dir/assemblies_pass/*.fa > $output_dir/assemblies_pass.fasta
    singularity exec $nc_d nextclade --input-fasta $output_dir/assemblies_pass.fasta --output-csv $output_dir/nextclade_report_clearlabs.csv
    mv params_clearlabs.yaml $output_dir

  else
    ###### run normal pipeline
    echo "run normal pipeline"
    nextflow run ~/Applications/daytona/flaq_sc2_humanclean2.nf -params-file params.yaml

    sort ./output/*/report.txt | uniq > ./output/sum_report.txt
    sed -i '/sampleID\treference/d' ./output/sum_report.txt
    sed -i '1i sampleID\treference\tstart\tend\tnum_raw_reads\tnum_clean_reads\tnum_mapped_reads\tpercent_mapped_clean_reads\tcov_bases_mapped\tpercent_genome_cov_map\tmean_depth\tmean_base_qual\tmean_map_qual\tassembly_length\tnumN\tpercent_ref_genome_cov\tVADR_flag\tQC_flag\tpangolin_version\tlineage\tSOTC' ./output/sum_report.txt

    cat ./output/assemblies/*.fa > ./output/assemblies.fasta
    singularity exec $nc_d nextclade --input-fasta ./output/assemblies.fasta --output-csv ./output/nextclade_report.csv
    mv params_clearlabs.yaml $output
  fi
  
else
  #---------------------------------------
  # Generating params_clearlabs.yaml file
  #---------------------------------------

  mkdir $output_dir

  get_abs_filename() {
    # $1 : relative filename
    echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
  }

  fastq="/data/fastq"
  alignment="/data/bam"
  fasta="/data/fasta"


  base=$(get_abs_filename $from_dir)
  outbase=$(get_abs_filename $output_dir)

  input="\"${base}${fastq}\""
  bams="\"${base}${alignment}\""
  assemblies="\"${base}${fasta}\""
  output="\"${outbase}\""
  echo "$output"
  sotc="\"$sotc\""

  pangolin_docker="\"$pangolin_docker\""
  pangolin_version="\"$pangolin_version\""
  pangolin_data_version="\"$pdata\""

  vadr_docker="\"$vadr_docker\""
  samtools_docker="\"$samtools_docker\""
  nc_d=$nextclade_docker
  nextclade_docker="\"$nextclade_docker\""

  touch params_clearlabs.yaml
  echo "# The parameters "input", "bams", "assemblies" and "output" are the absolute paths of the folders. Do not include the "/" at the end of the paths." >> params_clearlabs.yaml
  echo "input : "$input | tee -a params_clearlabs.yaml
  echo "bams : "$bams | tee -a params_clearlabs.yaml
  echo "assemblies : "$assemblies | tee -a params_clearlabs.yaml
  echo "output : "$output | tee -a params_clearlabs.yaml
  echo ""| tee -a params_clearlabs.yaml
  echo "# comma separated list of SOTCs to screen, default as "S:L452R,S:E484K""| tee -a params_clearlabs.yaml
  echo "sotc : "$sotc | tee -a params_clearlabs.yaml
  echo ""| tee -a params_clearlabs.yaml
  echo "#Do not change the parameter settings below:" | tee -a params_clearlabs.yaml
  echo "pangolin_docker : "$pangolin_docker | tee -a params_clearlabs.yaml
  echo "pangolin_version : "$pangolin_version | tee -a params_clearlabs.yaml
  echo "pangolin_data_version : "$pangolin_data_version | tee -a params_clearlabs.yaml

  echo ""| tee -a params_clearlabs.yaml
  echo "vadr_docker : "$vadr_docker | tee -a params_clearlabs.yaml
  echo "samtools_docker : "$samtools_docker | tee -a params_clearlabs.yaml
  echo "nextclade_docker : "$nextclade_docker | tee -a params_clearlabs.yaml

  #---------------------------
  # Running Daytona Pipeline
  #---------------------------

  if [ "$daytona_variation" = "clearlabs" ]; then
    ##### run clearlabs pipeline
    echo "run clearlabs pipeline"
    nextflow run ~/Applications/daytona/flaq_sc2_clearlabs2.nf -params-file params_clearlabs.yaml

    sort $output_dir/*/report.txt | uniq > $output_dir/sum_report.txt
    sed -i '/sampleID\treference/d' $output_dir/sum_report.txt
    sed -i '1i sampleID\treference\tstart\tend\tnum_clean_reads\tnum_mapped_reads\tpercent_mapped_clean_reads\tcov_bases_mapped\tpercent_genome_cov_map\tmean_depth\tmean_base_qual\tmean_map_qual\tassembly_length\tnumN\tpercent_ref_genome_cov\tVADR_flag\tQC_flag\tpangolin_version\tlineage\tSOTC' $output_dir/sum_report.txt
    
    awk 'BEGIN {FS=OFS="\t"} {print $1,$8,$10,$14,$15,$17,$18,$19}' $output_dir/sum_report.txt > $output_dir/mid_file
    paste -d'\t' $multi_mets $output_dir/mid_file | awk 'BEGIN {FS=OFS="\t"} {print $3,$4,$5,$6,$7,$8,$9,$10,$2}' > $output_dir/daytona_report.tsv
    rm $output_dir/mid_file

    cat $output_dir/assemblies_pass/*.fa > $output_dir/assemblies_pass.fasta
    singularity exec $nc_d nextclade --input-fasta $output_dir/assemblies_pass.fasta --output-csv $output_dir/nextclade_report_clearlabs.csv
    mv params_clearlabs.yaml $output_dir

  else
    ###### run normal pipeline
    echo "run normal pipeline"
    nextflow run ~/Applications/daytona/flaq_sc2_humanclean2.nf -params-file params.yaml

    sort ./output/*/report.txt | uniq > ./output/sum_report.txt
    sed -i '/sampleID\treference/d' ./output/sum_report.txt
    sed -i '1i sampleID\treference\tstart\tend\tnum_raw_reads\tnum_clean_reads\tnum_mapped_reads\tpercent_mapped_clean_reads\tcov_bases_mapped\tpercent_genome_cov_map\tmean_depth\tmean_base_qual\tmean_map_qual\tassembly_length\tnumN\tpercent_ref_genome_cov\tVADR_flag\tQC_flag\tpangolin_version\tlineage\tSOTC' ./output/sum_report.txt

    cat ./output/assemblies/*.fa > ./output/assemblies.fasta
    singularity exec $nc_d nextclade --input-fasta ./output/assemblies.fasta --output-csv ./output/nextclade_report.csv
    mv params_clearlabs.yaml $output
  fi
fi
