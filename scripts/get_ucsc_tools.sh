#!/usr/bin/env bash

mkdir ucsc_tools

TOOLS=(http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftUp, http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/axtChain http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chainMergeSort http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chainSplit http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chainNet http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/netChainSubset http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faSplit)

for i in ${TOOLS[@]}
do
    curl $i --output ucsc_tools/"${i##*/}"
done

chmod +x ucsc_tools/*

