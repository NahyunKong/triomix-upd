FROM elle72/tmp:vs2

WORKDIR /root

RUN R -e "install.packages('BiocManager')" \
    && R -e "BiocManager::install(c('DSS', 'BSSEQ'))"

RUN apt-get update && apt-get install -y \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    pkg-config \
    libfontconfig1-dev \
    && R -e "install.packages('tidyverse')"

RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    git \
    bc \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    python3 \
    python3-pip \
    && apt-get clean

RUN apt-get update && apt-get install -y \
    bedtools \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


RUN wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2 \
    && tar -xjf samtools-1.15.1.tar.bz2 \
    && cd samtools-1.15.1 \
    && make \
    && make install \
    && cd .. \
    && rm -rf samtools-1.15.1 samtools-1.15.1.tar.bz2


RUN git clone --recurse-submodules https://github.com/samtools/htslib.git && \
    cd htslib && \
    make && \
    make install




RUN rm -rf /var/lib/apt/lists/* htslib
