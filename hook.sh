. /anaconda/envs/azureml_py38/etc/conda/deactivate.d/libglib_deactivate.sh
. /anaconda/envs/azureml_py38/etc/conda/deactivate.d/tesseract_deactivate.sh
export PATH="/mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/hydromt_v1_somerset/.pixi/envs/compass-v1/bin:/home/azureuser/.pixi/bin:/home/azureuser/.vscode-server/data/User/globalStorage/github.copilot-chat/debugCommand:/home/azureuser/.vscode-server/data/User/globalStorage/github.copilot-chat/copilotCli:/home/azureuser/.vscode-server/bin/0958016b2af9f09bb4257e0df4a95e2f90590f9f/bin/remote-cli:/home/azureuser/bin:/home/azureuser/.local/bin:/home/azureuser/.pixi/bin:/home/azureuser/.vscode-server/bin/0958016b2af9f09bb4257e0df4a95e2f90590f9f/bin:/home/azureuser/bin:/home/azureuser/.local/bin:/home/azureuser/.juliaup/bin:/home/azureuser/.pixi/bin:/anaconda/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/home/azureuser/.dotnet/tools:/home/azureuser/.dotnet/tools:/home/azureuser/.vscode-server/extensions/ms-python.debugpy-2026.6.0-linux-x64/bundled/scripts/noConfigScripts"
export CONDA_SHLVL=7
export CONDA_ENV_SHLVL_7_CONDA_PREFIX=/anaconda/envs/azureml_py38
export CONDA_PREFIX=/mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/hydromt_v1_somerset/.pixi/envs/compass-v1
export PIXI_PROJECT_VERSION=0.1.0
export PIXI_IN_SHELL=1
export PIXI_PROJECT_ROOT=/mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/hydromt_v1_somerset
export PIXI_PROJECT_NAME=COMPASS_UC3
export PIXI_PROJECT_MANIFEST=/mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/hydromt_v1_somerset/pixi.toml
export PIXI_EXE=/home/azureuser/.pixi/bin/pixi
export CONDA_ENV_SHLVL_7_CONDA_DEFAULT_ENV=azureml_py38
export CONDA_DEFAULT_ENV=COMPASS_UC3:compass-v1
export PIXI_ENVIRONMENT_NAME=compass-v1
export PIXI_ENVIRONMENT_PLATFORMS='win-64,linux-64'
export PIXI_PROMPT='(COMPASS_UC3:compass-v1) '
. /mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/hydromt_v1_somerset/.pixi/envs/compass-v1/etc/conda/activate.d/gdal-activate.sh
. /mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/hydromt_v1_somerset/.pixi/envs/compass-v1/etc/conda/activate.d/libarrow_activate.sh
. /mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/hydromt_v1_somerset/.pixi/envs/compass-v1/etc/conda/activate.d/libglib_activate.sh
. /mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/hydromt_v1_somerset/.pixi/envs/compass-v1/etc/conda/activate.d/libxml2-split_activate.sh
. /mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/hydromt_v1_somerset/.pixi/envs/compass-v1/etc/conda/activate.d/proj4-activate.sh
export PYTHONNOUSERSITE=1
source /mnt/batch/tasks/shared/LS_root/mounts/clusters/ci-5-ukcrcompassflood/code/Users/eloise.matthews/hydromt_v1_somerset/.pixi/envs/compass-v1/share/bash-completion/completions/*

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

export PS1="(COMPASS_UC3:compass-v1) ${PS1:-}"
