#!/bin/bash

#EXIT IF INSTALLATION HAS WARNING OR ERROR
set -e
set -u
set -o pipefail

#OUTPUT TO LOG
exec >  >(tee -ia setup.log)
exec 2> >(tee -ia setup.log >&2)

#VARIABLES
TMPFILE=`mktemp`
PWD=`pwd`
now=`date`

#SOFTWARES
HISAT2_URL=ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
HISAT2_V=hisat2-2.1.0
BBMAP_URL=https://downloads.sourceforge.net/project/bbmap/BBMap_38.05.tar.gz
BBMAP_V=bbmap
LAST_URL=http://last.cbrc.jp/last-941.zip
LAST_V=last-941

#BAM BED PROCESSING
BEDTOOLS2_URL=https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz
BEDTOOLS2_V=bedtools2
SAMTOOLS_URL=https://downloads.sourceforge.net/project/samtools/samtools/1.8/samtools-1.8.tar.bz2
SAMTOOLS_V=samtools-1.8
SAMBAMBA_URL=https://github.com/biod/sambamba/releases/download/v0.6.8/sambamba-0.6.8-linux-static.gz
SAMBAMBA_V=sambamba

#UCSC
BG2BW_URL=http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bedGraphToBigWig
BG2BW_V=bedGraphToBigWig
FASR_URL=http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/faSomeRecords
FASR_V=faSomeRecords
FASZ_URL=http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/faSize
FASZ_V=faSize

#R PACKAGES
PLOTRIX_URL=https://cran.r-project.org/src/contrib/Archive/plotrix/plotrix_3.7-3.tar.gz
PLOTRIX_V=plotrix_3.7-3.tar.gz
RCB_URL=https://cran.r-project.org/src/contrib/Archive/RColorBrewer/RColorBrewer_1.0-5.tar.gz
RCB_V=RColorBrewer_1.0-5.tar.gz

#GENOME DATA
HG19_URL=https://www.dropbox.com/s/e1xwzye9hieewxz/human_g1k_v37.fasta.gz
HG19_V=hg19_g1k.fasta
MTRCRS_V=human_mt_rCRS.fasta
HG19S_V=hg19_g1k.size
M38_URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz
M38_V=mm10.fasta
MTM_V=mouse_mt.fasta
M38S_V=mm10.size

FOLDER=yes
HISAT=yes
BEDTOOLS=yes
BBMAP=yes
LAST=yes
SAMTOOL=yes
SAMBAMBA=yes
UCSC=yes
HUMAN=yes
HINDEX=yes
MOUSE=yes
MINDEX=yes
PLOTRIX=yes
RCB=yes
#echo -e "\x1b[97;41m Number of threads to use for building hisat index (your system info is given below): \x1b[m"
#lscpu | grep -E '^Thread|^Core|^Socket|^CPU\('
#read THREADS;

echo -e "\n\n\x1b[30;43m $now: Ok! based on your selection it may take upto two hours to complete all the download and compiling!\x1b[m"
#CLEAN DATA
if [ "$FOLDER" == "yes" ]; then
  echo -e "bam\nbin\nbw\ngenome\nindel\nlog\nplot\ntab" > tmp.list
  while read -r dir; do
    if [[ -d $dir ]]; then
      rm -r $dir
      mkdir $dir
    else
      mkdir $dir
    fi  
  done < tmp.list
  rm tmp.list
fi

#DOWNLOAD AND CONFIGURE SOFTWARES REQUIRED
echo -e "\x1b[97;41m $now:Installing Hisat2: \x1b[m"
if [ "$HISAT" == "yes" ]; then
wget -O $TMPFILE $HISAT2_URL
unzip -d $PWD $TMPFILE
rm $TMPFILE
mv $HISAT2_V bin/hisat2
fi

echo -e "\x1b[97;41m $now:Installing Bedtools: \x1b[m"
if [ "$BEDTOOLS" == "yes" ]; then
wget -O $TMPFILE $BEDTOOLS2_URL
tar xfvz $TMPFILE -C $PWD
rm $TMPFILE
make -C $BEDTOOLS2_V/
mv $BEDTOOLS2_V bin/bedtools2
fi

echo -e "\x1b[97;41m $now:Installing BBMap: \x1b[m"
if [ "$BBMAP" == "yes" ]; then
wget -O $TMPFILE $BBMAP_URL
tar xfvz $TMPFILE -C $PWD/bin/
rm $TMPFILE
fi


echo -e "\x1b[97;41m $now:Installing LAST: \x1b[m"
if [ "$LAST" == "yes" ]; then
wget -O $TMPFILE $LAST_URL
unzip -d $PWD $TMPFILE
rm $TMPFILE
make -C $LAST_V/src/
mv $LAST_V bin/last
fi

echo -e "\x1b[97;41m $now:Installing Samtools: \x1b[m"
if [ "$SAMTOOL" == "yes" ]; then
wget -O $TMPFILE $SAMTOOLS_URL
tar xfvj $TMPFILE -C $PWD/
rm $TMPFILE
cd $SAMTOOLS_V
./configure
make
cd ..
mv $SAMTOOLS_V bin/samtools
fi

echo -e "\x1b[97;41m $now:Installing Sambamba: \x1b[m"
if [ "$SAMBAMBA" == "yes" ]; then
wget -O $TMPFILE $SAMBAMBA_URL
gunzip -c $TMPFILE > bin/$SAMBAMBA_V
chmod 755 bin/$SAMBAMBA_V
rm $TMPFILE
fi

#GET UCSC UTILITIES
echo -e "\x1b[97;41m $now:Installing UCSC: \x1b[m"
if [ "$UCSC" == "yes" ]; then
wget -O $TMPFILE $BG2BW_URL
mv $TMPFILE bin/$BG2BW_V
chmod 755 bin/$BG2BW_V

wget -O $TMPFILE $FASR_URL
mv $TMPFILE bin/$FASR_V
chmod 755 bin/$FASR_V

wget -O $TMPFILE $FASZ_URL
mv $TMPFILE bin/$FASZ_V
chmod 755 bin/$FASZ_V
fi

#DOWNLOAD HUMAN GENOME AND EXTRACT MT GENOME
echo -e "\x1b[97;41m $now:Download and index human genome: \x1b[m"
if [ "$HUMAN" == "yes" ]; then
  wget -O $TMPFILE $HG19_URL
  gunzip -c $TMPFILE > genome/$HG19_V
  bin/$FASZ_V -detailed genome/$HG19_V|egrep -v 'GL|MT' > genome/$HG19S_V
  rm $TMPFILE

  echo 'MT' > tmp.txt
  bin/$FASR_V genome/$HG19_V tmp.txt genome/$MTRCRS_V
  rm tmp.txt
fi

#BUILD HUMAN INDEXES
if [ "$HINDEX" == "yes" ]; then
  #bin/hisat2/hisat2-build -p $THREADS genome/$HG19_V genome/hg19_g1k
  bin/last/src/lastdb -uNEAR  genome/human_mt_rCRS genome/$MTRCRS_V
  bin/samtools/samtools faidx genome/$HG19_V
  bin/samtools/samtools faidx genome/$MTRCRS_V
fi

#DOWNLOAD MOUSE GENOME AND EXTRACT MT GENOME
echo -e "\x1b[97;41m $now:Download and index mouse genome: \x1b[m"
if [ "$MOUSE" == "yes" ]; then
  wget -O $TMPFILE $M38_URL
  gunzip -c $TMPFILE > genome/$M38_V
  bin/$FASZ_V -detailed genome/$M38_V|fgrep NC|fgrep -v 'NC_005089.1' > genome/$M38S_V
  rm $TMPFILE

  echo 'NC_005089.1' > tmp.txt
  bin/$FASR_V genome/$M38_V tmp.txt genome/$MTM_V
  rm tmp.txt
fi

#BUILD MOUSE INDEXES
if [ "$MINDEX" == "yes" ]; then
  #bin/hisat2/hisat2-build -p $THREADS genome/$M38_V genome/mm10
  bin/last/src/lastdb -uNEAR  genome/mouse_mt genome/$MTM_V
  bin/samtools/samtools faidx genome/$M38_V
  bin/samtools/samtools faidx genome/$MTM_V
fi

#INSTALL PLOTRIX
echo -e "\x1b[97;41m $now:Install R packages: \x1b[m"
if [ "$PLOTRIX" == "yes" ]; then
  wget $PLOTRIX_URL
  R CMD INSTALL -l bin $PLOTRIX_V
  rm $PLOTRIX_V
fi

#INSTALL RCOLORBREWER
if [ "$RCB" == "yes" ]; then
  wget $RCB_URL
  R CMD INSTALL -l bin $RCB_V
  rm $RCB_V
fi

echo -e "\n\n\x1b[30;43m $now: Setup finished!\x1b[m"
