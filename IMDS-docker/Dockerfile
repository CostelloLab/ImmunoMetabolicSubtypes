# Full contents of Dockerfile

FROM continuumio/miniconda3
LABEL description="Base docker image with conda and util libraries"
ARG ENV_NAME="IM-DS"

# Install the conda environment
RUN conda create -n ${ENV_NAME} -y


# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/${ENV_NAME}/bin:$PATH
