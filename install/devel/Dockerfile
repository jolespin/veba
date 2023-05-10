FROM ubuntu:latest

ARG ENV_NAME

SHELL ["/bin/bash", "-c"]

WORKDIR /root/

# Install required packages
RUN apt-get update && \
    apt-get install -y coreutils wget bzip2 ca-certificates curl git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Retrieve VEBA repository
RUN mkdir -p veba/
COPY ./install/ veba/install/
COPY ./src/ veba/src/
COPY ./VERSION veba/VERSION
COPY ./LICENSE veba/LICENSE

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    /opt/conda/bin/conda init bash && \
    /opt/conda/bin/conda config --add channels jolespin && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    /opt/conda/bin/conda install -c conda-forge mamba && \
    /opt/conda/bin/mamba init bash && \
    /opt/conda/bin/conda clean -afy

# use this file structure for. installing miniconda
ENV PATH /opt/conda/bin:$PATH

# Create environment
RUN mamba env create -n ${ENV_NAME} -f veba/install/environments/${ENV_NAME}.yml

# Add environment scripts to 
RUN bash veba/install/update_environment_scripts.sh veba/

# Set up environment
# ENTRYPOINT ["bash", "-c", "/opt/conda/bin/conda activate ${ENV_NAME}"]