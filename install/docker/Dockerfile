# v2023.6.6
# =================================
# Miniconda3
# =================================

FROM continuumio/miniconda3

ARG ENV_NAME

SHELL ["/bin/bash","-l", "-c"]

WORKDIR /root/

# Data
RUN mkdir -p /volumes/input
RUN mkdir -p /volumes/output
RUN mkdir -p /volumes/database


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

