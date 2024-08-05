FROM ubuntu:20.04

MAINTAINER Rick Wertenbroek <rick.wertenbroek@unil.ch>

ENV DEBIAN_FRONTEND noninteractive

# Install required software and clean as not to make the layer dirty
RUN apt-get update && apt-get -y upgrade && apt-get install -y \
	apt-utils curl gnupg gcc g++ make autoconf && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update && apt-get -y upgrade && apt-get install -y \
	git zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update && apt-get -y upgrade && apt-get install -y \
	bcftools samtools unzip && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# https://github.com/samtools/htslib/blob/f68c84ff46a63fde05ae2bdf3ebdec568c9c4613/cram/cram_io.c#L2923
RUN curl -LJO https://github.com/rwk-unil/sapphire/releases/download/v1.0.0-paper/hts-ref.zip && \
	mkdir -p /root/.cache && mv hts-ref.zip /root/.cache/ && cd /root/.cache/ && unzip hts-ref.zip && rm -f hts-ref.zip

# Clone Source Code Repository
RUN mkdir -p /usr/src/ && \
    cd /usr/src/ && \
    git clone https://github.com/rwk-unil/sapphire.git && \
    cd /usr/src/sapphire && \
    cd Docker && \
    chmod +x *.sh && \
    ./prepare.sh

RUN cd /usr/src/sapphire/Docker && \
    ./build_install.sh

# Work in this temporary directory
WORKDIR /tmp/work

CMD echo "Run with the following command : docker run <tag> <command> [args]"