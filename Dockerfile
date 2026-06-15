# Dockerfile for ModelPipelineIO CWL workflows
# Uses micromamba for faster dependency resolution and smaller image size

FROM mambaorg/micromamba:latest

# Create and activate the modelpio environment
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean -a --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install the ModelPipelineIO package in development mode
WORKDIR /app
COPY --chown=$MAMBA_USER . /app

# RUN chown -R $MAMBA_USER:$MAMBA_USER /app
RUN pip install -e .

# Set working directory to CWL scripts
WORKDIR /app/cwl

# Set the path to include the micromamba directory
ENV PATH=/opt/mamba/bin:$PATH

# Default command (can be overridden)
CMD [\"cwltool\", \"--help\"]
