# v2023.8.22
# =================================
FROM mambaorg/micromamba:1.4.9

ARG ENV_NAME

SHELL ["/usr/local/bin/_dockerfile_shell.sh"]

WORKDIR /tmp/

# Data
USER root
RUN mkdir -p /volumes/
RUN mkdir -p /volumes/input
RUN mkdir -p /volumes/output
RUN mkdir -p /volumes/database

# Retrieve VEBA repository
RUN mkdir -p veba/
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER ./install/ veba/install/
COPY --chown=$MAMBA_USER:$MAMBA_USER ./src/ veba/src/
COPY --chown=$MAMBA_USER:$MAMBA_USER ./VERSION veba/VERSION
COPY --chown=$MAMBA_USER:$MAMBA_USER ./LICENSE veba/LICENSE

# Install dependencies
RUN micromamba install -y -n base -f veba/install/environments/${ENV_NAME}.yml && \ 
    micromamba clean -a -y -f

# Add environment scripts to environment bin
RUN cp -rf veba/src/* /opt/conda/bin/ && \
    ln -sf /opt/conda/bin/scripts/*.py /opt/conda/bin/ && \
    ln -sf /opt/conda/bin/scripts/*.r /opt/conda/bin/


ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
