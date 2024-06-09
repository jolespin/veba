# v2024.4.29
# =================================
FROM --platform=linux/amd64 mambaorg/micromamba:latest

ARG ENV_NAME

SHELL ["/usr/local/bin/_dockerfile_shell.sh"]

WORKDIR /tmp/

# Data
USER root
RUN mkdir -p /volumes/
RUN mkdir -p /volumes/input
RUN mkdir -p /volumes/output
RUN mkdir -p /volumes/database
RUN mkdir -p /volumes/workspace

# Retrieve VEBA repository
RUN mkdir -p veba/
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER ./install/ veba/install/
COPY --chown=$MAMBA_USER:$MAMBA_USER ./bin/ veba/bin/
COPY --chown=$MAMBA_USER:$MAMBA_USER ./VERSION veba/VERSION
COPY --chown=$MAMBA_USER:$MAMBA_USER ./LICENSE veba/LICENSE
RUN chmod -R 777 veba/bin/

# Fix: 
# critical libmamba Download error (28) Timeout was reached
# [https://conda.anaconda.org/bioconda/linux-64/blast-2.13.0-hf3cf87c_0.tar.bz2]
# Operation too slow. Less than 30 bytes/sec transferred the last 60 seconds

ARG MAMBA_NO_LOW_SPEED_LIMIT=1

# Build .condarc
# RUN echo "channel_priority: flexible" >> ~/.condarc
# RUN echo "channels:" >> ~/.condarc
# RUN echo "  - conda-forge" >> ~/.condarc
# RUN echo "  - bioconda" >> ~/.condarc
# RUN echo "  - jolespin" >> ~/.condarc

# Install dependencies
RUN micromamba install -y -n base -f veba/install/environments/${ENV_NAME}.yml && \ 
    micromamba clean -a -y -f

# Add environment scripts to environment bin
RUN cp -rf veba/bin/* /opt/conda/bin/ && \
    ln -sf /opt/conda/bin/scripts/*.py /opt/conda/bin/ && \
    ln -sf /opt/conda/bin/scripts/*.r /opt/conda/bin/

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
