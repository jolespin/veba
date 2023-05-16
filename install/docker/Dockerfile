# v2023.5.11
# # =================================
# # Ubuntu
# # =================================
# FROM ubuntu:latest

# ARG ENV_NAME

# SHELL ["/bin/bash", "-c"]

# WORKDIR /root/

# # Install required packages
# RUN apt-get update && \
#     apt-get install -y coreutils wget bzip2 ca-certificates curl git && \
#     apt-get clean && \
#     rm -rf /var/lib/apt/lists/*

# # Retrieve VEBA repository
# RUN mkdir -p veba/
# COPY ./install/ veba/install/
# COPY ./src/ veba/src/
# COPY ./VERSION veba/VERSION
# COPY ./LICENSE veba/LICENSE

# # Install Miniconda
# RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
#     bash miniconda.sh -b -p /opt/conda && \
#     rm ~/miniconda.sh && \
#     ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
#     /opt/conda/bin/conda init bash && \
#     /opt/conda/bin/conda config --add channels jolespin && \
#     /opt/conda/bin/conda config --add channels bioconda && \
#     /opt/conda/bin/conda config --add channels conda-forge && \
#     /opt/conda/bin/conda update conda -y && \
#     /opt/conda/bin/conda install -c conda-forge mamba -y && \
#     /opt/conda/bin/mamba init bash && \
#     /opt/conda/bin/conda clean -afy

# =================================
# Miniconda3
# =================================

FROM continuumio/miniconda3

ARG ENV_NAME

SHELL ["/bin/bash","-l", "-c"]

WORKDIR /root/

# Data
RUN mkdir -p /data/

# Retrieve VEBA repository
RUN mkdir -p veba/
COPY ./install/ veba/install/
COPY ./src/ veba/src/
COPY ./VERSION veba/VERSION
COPY ./LICENSE veba/LICENSE

# Install Miniconda
RUN /opt/conda/bin/conda init bash && \
    /opt/conda/bin/conda config --add channels jolespin && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    /opt/conda/bin/conda update conda -y && \
    # /opt/conda/bin/conda install -c conda-forge mamba -y && \ # Mamba adds about 450 MB to image
    # /opt/conda/bin/mamba init bash && \
    /opt/conda/bin/conda clean -afy

# =================================

# Add conda bin to path
ENV PATH /opt/conda/bin:$PATH

# Create environment
RUN conda env create -n ${ENV_NAME} -f veba/install/environments/${ENV_NAME}.yml
# RUN mamba env create -n ${ENV_NAME} -f veba/install/environments/${ENV_NAME}.yml

# Add environment scripts to environment bin
RUN /bin/bash veba/install/update_environment_scripts.sh veba/

# # Add contents to path
# ENV PATH /opt/conda/envs/${ENV_NAME}/bin:$PATH

# Set up environment
RUN echo "conda activate ${ENV_NAME}" >> ~/.bashrc

# Set entrypoint to bash
ENTRYPOINT ["bash", "-l", "-c"]

