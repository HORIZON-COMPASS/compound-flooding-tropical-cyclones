# Image for the COMPASS Snakemake workflow on AzureML AmlCompute.
#
# The conda spec alone is not enough: the workflow shells out to two native solvers.
#   - Wflow  -> Julia + Wflow.jl. The bundled wflow_cli.exe is a Windows PE binary and
#               cannot run on Linux, so snakefile_wflow.smk always falls through to julia.
#   - SFINCS -> apptainer running 02_Models/00_executables/sfincs-cpu_v2.4.1.sif.
#
# Build context is 00_execution_script_examples/aml_image/, holding only this file and
# aml_env_compass-v1.yml. Keep anything else out: editing any file in the context triggers
# a full image rebuild on the next submit.

FROM mcr.microsoft.com/azureml/openmpi4.1.0-ubuntu22.04:latest

# --------------------------------------------------------------------------------------
# Apptainer, to run the SFINCS .sif that already lives on the file share.
# tzdata provides /etc/localtime, which apptainer bind-mounts into every container by
# default; the base image lacks it and `apptainer run` fails with "mount source
# /etc/localtime doesn't exist".
# --------------------------------------------------------------------------------------
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      software-properties-common curl ca-certificates tzdata && \
    add-apt-repository -y ppa:apptainer/ppa && \
    apt-get update && \
    apt-get install -y --no-install-recommends apptainer && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# --------------------------------------------------------------------------------------
# Conda environment (pinned to .pixi/envs/compass-v1 — see aml_env_compass-v1.yml).
# --------------------------------------------------------------------------------------
COPY aml_env_compass-v1.yml /tmp/aml_env_compass-v1.yml
RUN conda env create -p /azureml-envs/compass-v1 -f /tmp/aml_env_compass-v1.yml && \
    rm -rf "$HOME/.cache/pip" && \
    conda clean -aqy

ENV PATH=/azureml-envs/compass-v1/bin:$PATH \
    CONDA_PREFIX=/azureml-envs/compass-v1 \
    LD_LIBRARY_PATH=/azureml-envs/compass-v1/lib:$LD_LIBRARY_PATH

# --------------------------------------------------------------------------------------
# Julia + Wflow.jl, matching the compute instance (julia 1.12.4, Wflow 1.0.4).
# snakefile_wflow.smk calls `julia --project=~/.julia/environments/wflow1`; julia expands
# the tilde, so the job command symlinks $HOME/.julia/environments/wflow1 to this project.
# --------------------------------------------------------------------------------------
ARG JULIA_VERSION=1.12.4
ARG JULIA_SERIES=1.12
ARG WFLOW_VERSION=1.0.4

ENV JULIA_DEPOT_PATH=/opt/julia-depot
ENV WFLOW_PROJECT=/opt/julia-depot/environments/wflow1

# `julia` is a wrapper, not a symlink: the conda lib dir on LD_LIBRARY_PATH shadows Julia's
# bundled libcurl/libstdc++ etc. and Julia segfaults (exit 139) as soon as it loads Pkg.
RUN curl -fsSL "https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_SERIES}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz" \
      | tar -xz -C /opt && \
    printf '#!/bin/sh\nexec env -u LD_LIBRARY_PATH /opt/julia-%s/bin/julia "$@"\n' "${JULIA_VERSION}" \
      > /usr/local/bin/julia && \
    chmod 755 /usr/local/bin/julia && \
    julia --version

# Precompile at build time so job startup does not pay for it.
RUN mkdir -p "${WFLOW_PROJECT}" && \
    julia --project="${WFLOW_PROJECT}" -e \
      "using Pkg; Pkg.add(name=\"Wflow\", version=\"${WFLOW_VERSION}\"); Pkg.precompile(); using Wflow" && \
    chmod -R a+rX /opt/julia-depot

ENV HDF5_USE_FILE_LOCKING=FALSE
