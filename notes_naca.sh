
# UCSC genome browser exon locations to determine what mapping was used with R
[1] "chr12:57106211-57111211" "chr12:57106570-57111570" "chr12:57106846-57111846" "chr12:57107321-57112321"
[5] "chr12:57108146-57113146" "chr12:57108386-57113386" "chr12:57109655-57114655" "chr12:57118236-57123236"
[9] "chr12:57118746-57123746"
# turns out exon 2 ensembl mapping was used. didn't try to differentiate for missing domain, so actually most reads were kept

#naca_investigation.R: 1-00425, 1-01026, 1-05824
# IGV bam URLs obtained from url.hlhs.bam.list
cd /Users/frichter/Documents/phd/variancepartitionpcgc/analysisExon/DE_Exon_F30/L2HLHS 
grep "00425\|01026\|05824" url.hlhs.bam.list > url.hlhs.SPLICED.bam.list
# IGV: File → load from file
# obtained by ctrl + F sapien at ensembl distributed annotation server:
# http://das.ensembl.org/das/sources
# http://das.ensembl.org/das/Homo_sapiens.GRCh37.transcript/features
# IGV: File → load from DAS

# when screenshotting Sashimi plot be sure to capture name

grep -v "00425\|01026\|05824" url.hlhs.bam.list > url.hlhs.NOTSPLICED.bam.list
head -n25 url.hlhs.not.bam.list > url.hlhs.not.5.bam.list
tail -n25 url.hlhs.not.bam.list > url.hlhs.not.6.bam.list
# url.hlhs.not.5.bam.list corrupted? Memory too low