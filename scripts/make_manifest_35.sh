# Makes manifest file for input into qiime 2 (2020.8)
# specifically for MiSeqV1V3_35 run (qi's demultiplexed reads)

pathdemux="/Users/amycampbell/Documents/IowaWoundData2021/MiSeqV1V3_35/demultiplexed/MiSeqV1V3_35_barcode_"
extfwd="_1.fastq.gz"
extrev="_2.fastq.gz"
manifestpath="/Users/amycampbell/Desktop/GriceLabGit/IowaWound/mappings/Manifest35.tsv"
mappingfile="/Users/amycampbell/Desktop/GriceLabGit/IowaWound/mappings/IA_woundpain_mapping_35_2021.csv"
intermediates="/Users/amycampbell/Desktop/GriceLabGit/IowaWound/intermediates/"

# Make folder for intermediate files
####################################
mkdir -p $intermediates

# make header
echo "sample-id	forward-absolute-filepath	reverse-absolute-filepath" >  "${intermediates}header.txt"

# make list of sample IDs 
tail -n +2 $mappingfile | cut -d, -f1 > "${intermediates}samplelist.txt"

# make fwd list of files 
awk -v pathvarawk="$pathdemux" 'BEGIN {FS="\t";}; {print pathvarawk $0 "_1.fastq.gz"}' "${intermediates}samplelist.txt" > "${intermediates}samplelist_fwd.txt"

# make rev list of files
awk -v pathvarawk="$pathdemux" 'BEGIN {FS="\t";}; {print pathvarawk $0 "_2.fastq.gz"}' "${intermediates}samplelist.txt" > "${intermediates}samplelist_rev.txt"


# join horizontally
#####################
paste "${intermediates}samplelist.txt" "${intermediates}samplelist_fwd.txt" > "${intermediates}samplelist_missingrev.txt"

paste "${intermediates}samplelist_missingrev.txt" "${intermediates}samplelist_rev.txt" > "${intermediates}samplelist_missingheader.txt"

# add header
############
cat "${intermediates}header.txt" "${intermediates}samplelist_missingheader.txt" > $manifestpath
