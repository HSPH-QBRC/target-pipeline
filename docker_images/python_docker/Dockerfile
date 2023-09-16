FROM ubuntu:jammy

SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y \
    build-essential \
    git \
    curl \
    python3 \
    python3-pip


# Install all software under /opt/software:
RUN mkdir -p /opt/software 
RUN mkdir -p /opt/software/resources
RUN mkdir -p /opt/software/scripts

# Add the various python requirements and script files
ADD requirements.txt /opt/software/

# for mapping ENSG to symbols
RUN curl ftp://ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip -o /opt/software/resources/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip

# Install python deps
RUN pip3 install -r /opt/software/requirements.txt

ADD map_ensg_to_symbol.py /opt/software/scripts/
ADD cohort_selection.py /opt/software/scripts/