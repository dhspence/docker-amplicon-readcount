FROM registry.gsc.wustl.edu/genome/gatk-3.6:1
MAINTAINER David Spencer <dspencer@wustl.edu>

LABEL \
    description="GATK plus custom python script to tag HaloplexHS VCF file with amplicon information"

COPY addAmpliconInfoAndCountReads.py /usr/local/bin/
