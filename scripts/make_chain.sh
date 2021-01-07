#!/usr/bin/env bash

## Download the following: (or run the get_ucsc_tools.sh script)

# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/axtChain
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chainMergeSort
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chainSplit
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chainNet
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/netChainSubset
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faSplit

## Required samtools in $PATH

## Steps adapted from: http://genomewiki.ucsc.edu/index.php/LiftOver_Howto

rm -rf net
rm -rf tmp
rm -rf chain

mkdir -p tmp/2bit
mkdir -p tmp/psl
mkdir -p tmp/gene_split
mkdir -p tmp/chain
mkdir -p tmp/build

CHROMOSOME_REF=$1
GENE_REF=$2



CHROMOSOME_BASE="${CHROMOSOME_REF##*/}"
GENE_BASE="${GENE_REF##*/}"

samtools faidx $CHROMOSOME_REF
cut -f1,2 $CHROMOSOME_REF.fai > tmp/build/$CHROMOSOME_BASE.genome

samtools faidx $GENE_REF
cut -f1,2 $GENE_REF.fai > tmp/build/$GENE_BASE.genome

ucsc_tools/faSplit size $GENE_REF 5000 tmp/gene_split/$GENE_BASE -lift=tmp/gene_split/lift

# Convert fasta to 2bit

ucsc_tools/faToTwoBit $CHROMOSOME_REF tmp/2bit/$CHROMOSOME_BASE.2bit
ucsc_tools/faToTwoBit $GENE_REF tmp/2bit/$GENE_BASE.2bit

# BLAT align

for fa in tmp/gene_split/*.fa
do
    ucsc_tools/blat tmp/2bit/$CHROMOSOME_BASE.2bit $fa tmp/psl/OLD."${fa##*/}".psl -tileSize=12 -minScore=100 -minIdentity=98 -fastMap
    ucsc_tools/liftUp -pslQ tmp/psl/"${fa##*/}".psl tmp/gene_split/lift warn tmp/psl/OLD."${fa##*/}".psl
    ucsc_tools/axtChain -linearGap=medium -psl tmp/psl/"${fa##*/}".psl tmp/2bit/$CHROMOSOME_BASE.2bit tmp/2bit/$GENE_BASE.2bit tmp/chain/"${fa##*/}".chain
done

ucsc_tools/chainMergeSort tmp/chain/*.chain | ucsc_tools/chainSplit chain stdin

mkdir net

cd chain

../ucsc_tools/chainNet *.chain ../tmp/build/$CHROMOSOME_BASE.genome ../tmp/build/$GENE_BASE.genome  ../net/final.net /dev/null

../ucsc_tools/netChainSubset ../net/final.net *.chain ../$CHROMOSOME_BASE-$GENE_BASE.chain

cd ..

rm -rf net
rm -rf tmp
rm -rf chain
