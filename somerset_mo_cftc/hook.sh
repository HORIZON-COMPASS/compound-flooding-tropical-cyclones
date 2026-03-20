. /anaconda/envs/azureml_py38/etc/conda/deactivate.d/libglib_deactivate.sh
. /anaconda/envs/azureml_py38/etc/conda/deactivate.d/tesseract_deactivate.sh
export PATH="/mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/cftc_fresh/somerset_mo_cftc/.pixi/envs/compass-wflow/bin:/home/azureuser/.juliaup/bin:/home/azureuser/.pixi/bin:/home/azureuser/.pixi/bin:/home/azureuser/.vscode-server/data/User/globalStorage/github.copilot-chat/debugCommand:/home/azureuser/.vscode-server/data/User/globalStorage/github.copilot-chat/copilotCli:/home/azureuser/.vscode-server/bin/585eba7c0c34fd6b30faac7c62a42050bfbc0086/bin/remote-cli:/home/azureuser/bin:/home/azureuser/.local/bin:/home/azureuser/.vscode-server/bin/585eba7c0c34fd6b30faac7c62a42050bfbc0086/bin:/home/azureuser/bin:/home/azureuser/.local/bin:/anaconda/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/home/azureuser/.dotnet/tools:/home/azureuser/.dotnet/tools:/home/azureuser/.vscode-server/extensions/ms-python.debugpy-2025.18.0-linux-x64/bundled/scripts/noConfigScripts"
export CONDA_SHLVL=11
export CONDA_ENV_SHLVL_11_CONDA_PREFIX=/anaconda/envs/azureml_py38
export CONDA_PREFIX=/mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/cftc_fresh/somerset_mo_cftc/.pixi/envs/compass-wflow
export CONDA_ENV_SHLVL_11_PIXI_PROJECT_ROOT=/mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/cftc_fresh/somerset_mo_cftc
export PIXI_PROJECT_ROOT=/mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/cftc_fresh/somerset_mo_cftc
export CONDA_ENV_SHLVL_11_PIXI_PROJECT_MANIFEST=/mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/cftc_fresh/somerset_mo_cftc/pixi.toml
export PIXI_PROJECT_MANIFEST=/mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/cftc_fresh/somerset_mo_cftc/pixi.toml
export CONDA_ENV_SHLVL_11_PIXI_PROJECT_NAME=COMPASS_UC3
export PIXI_PROJECT_NAME=COMPASS_UC3
export CONDA_ENV_SHLVL_11_PIXI_PROJECT_VERSION=0.1.0
export PIXI_PROJECT_VERSION=0.1.0
export CONDA_ENV_SHLVL_11_PIXI_IN_SHELL=1
export PIXI_IN_SHELL=1
export CONDA_ENV_SHLVL_11_PIXI_EXE=/home/azureuser/.pixi/bin/pixi
export PIXI_EXE=/home/azureuser/.pixi/bin/pixi
export CONDA_ENV_SHLVL_11_CONDA_DEFAULT_ENV=azureml_py38
export CONDA_DEFAULT_ENV=COMPASS_UC3:compass-wflow
export CONDA_ENV_SHLVL_11_PIXI_ENVIRONMENT_NAME=default
export PIXI_ENVIRONMENT_NAME=compass-wflow
export CONDA_ENV_SHLVL_11_PIXI_ENVIRONMENT_PLATFORMS='linux-64,win-64'
export PIXI_ENVIRONMENT_PLATFORMS='linux-64,win-64'
export CONDA_ENV_SHLVL_11_PIXI_PROMPT='(COMPASS_UC3) '
export PIXI_PROMPT='(COMPASS_UC3:compass-wflow) '
. /mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/cftc_fresh/somerset_mo_cftc/.pixi/envs/compass-wflow/etc/conda/activate.d/gdal-activate.sh
. /mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/cftc_fresh/somerset_mo_cftc/.pixi/envs/compass-wflow/etc/conda/activate.d/geotiff-activate.sh
. /mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/cftc_fresh/somerset_mo_cftc/.pixi/envs/compass-wflow/etc/conda/activate.d/libarrow_activate.sh
. /mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/cftc_fresh/somerset_mo_cftc/.pixi/envs/compass-wflow/etc/conda/activate.d/libglib_activate.sh
. /mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/cftc_fresh/somerset_mo_cftc/.pixi/envs/compass-wflow/etc/conda/activate.d/libxml2_activate.sh
. /mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/cftc_fresh/somerset_mo_cftc/.pixi/envs/compass-wflow/etc/conda/activate.d/proj4-activate.sh
. /mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/cftc_fresh/somerset_mo_cftc/.pixi/envs/compass-wflow/etc/conda/activate.d/udunits2-activate.sh
source /mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/cftc_fresh/somerset_mo_cftc/.pixi/envs/compass-wflow/share/bash-completion/completions/*

# shellcheck shell=bash
pixi() {
    local first_arg="${1-}"

    "${PIXI_EXE-}" "$@" || return $?

    case "${first_arg-}" in
    add | a | remove | rm | install | i)
        eval "$("$PIXI_EXE" shell-hook --change-ps1 false)"
        hash -r
        ;;
    esac || :

    return 0
}

export PS1="(COMPASS_UC3:compass-wflow) ${PS1:-}"
