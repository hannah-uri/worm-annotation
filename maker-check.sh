#!/bin/sh
#tests maker on wormbase data

mkdir ref
wget -P ref/ "ftp://ftp.wormbase.org/pub/wormbase/releases/WS274/species/c_elegans/PRJEB28388/c_elegans.PRJEB28388.WS274.genomic.fa.gz"
wget -P ref/ "ftp://ftp.wormbase.org/pub/wormbase/releases/WS274/species/c_elegans/PRJEB28388/c_elegans.PRJEB28388.WS274.mRNA_transcripts.fa.gz"
gunzip ref/c_elegans.PRJEB28388.WS274.genomic.fa.gz
gunzip ref/c_elegans.PRJEB28388.WS274.mRNA_transcripts.fa.gz

maker -CTL
cat maker_opts.ctl | sed 's/^genome=.*$/genome=ref\/c_elegans.PRJEB28388.WS274.genomic.fa /' | sed 's/^est=.*$/est=ref\/c_elegans.PRJEB28388.WS274.mRNA_transcripts.fa /' | sed 's/^model_org=.*$/model_org= /' > maker_opts.ctl.2
mv maker_opts.ctl.2 maker_opts.ctl

maker
