FROM ubuntu:20.04

MAINTAINER Rick Wertenbroek <rick.wertenbroek@unil.ch>

ENV DEBIAN_FRONTEND noninteractive

# Install required software and clean as not to make the layer dirty
RUN apt-get update && apt-get -y upgrade && apt-get install -y \
	apt-utils curl gnupg gcc g++ make autoconf && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update && apt-get -y upgrade && apt-get install -y \
	git zlib1g-dev libbz2-dev liblzma-dev && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update && apt-get -y upgrade && apt-get install -y \
	bcftools && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Clone Source Code Repository (with read-only token)
RUN mkdir -p /usr/src/ && \
    cd /usr/src/ && \
    git clone https://rwk-unil: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@github.com/rwk-unil/pp.git && \
    cd /usr/src/pp && \
    git submodule update --init --recursive xSqueezeIt && \
    cd xSqueezeIt && \
    cd htslib && \
    autoheader && \
    autoconf && \
    automake --add-missing; \
    ./configure && \
    make && \
    make install && \
    ldconfig && \
    cd .. && \
    git clone https://github.com/facebook/zstd.git && \
    cd zstd && \
    make && \
    cd .. && \
    cd .. && \
    cd pp_extractor && \
    make && \
    chmod +x pp_extract && \
    cp pp_extract /usr/local/bin/ && \
    chmod +x scripts/run_pp_extract.sh && \
    cp scripts/run_pp_extract.sh /usr/local/bin/run_pp_extract && \
    cd .. && \
    cd phase_caller && \
    make && \
    chmod +x phase_caller && \
    cp phase_caller /usr/local/bin/ && \
    cd /usr/src/ && \
    rm -rf pp

# Work in this temporary directory
WORKDIR /tmp/work

CMD echo "Run with the following command : docker run <tag> run_pp_extract [args]"