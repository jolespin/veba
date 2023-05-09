FROM ubuntu:latest

ARG ENV_NAME

#RUN chsh -s /bin/bash
SHELL ["/bin/bash", "-c"]

WORKDIR /root/

# Install required packages
RUN apt-get update && \
    apt-get install -y wget bzip2 ca-certificates curl git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    /opt/conda/bin/conda init bash && \
    /opt/conda/bin/conda clean -afy

# use this file structure for. installing miniconda
ENV PATH /opt/conda/bin:$PATH

# Install packages
COPY ./environments/${ENV_NAME}.yml ./environment.yml
COPY ./update_environment_scripts.sh ./

RUN conda env create -n ${ENV_NAME} -f environment.yml
RUN bash update_environment_scripts.sh ../