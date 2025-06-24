. "/u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/deactivate.d/gdal-deactivate.sh"
. "/u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/deactivate.d/geotiff-deactivate.sh"
. "/u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/deactivate.d/libglib_deactivate.sh"
. "/u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/deactivate.d/libxml2_deactivate.sh"
. "/u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/deactivate.d/proj4-deactivate.sh"
export PATH="/u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/bin:/opt/apps/juliaup/bin:/opt/apps/pixi:/u/morenodu/.pixi/bin:/u/morenodu/.vscode-server/cli/servers/Stable-18e3a1ec544e6907be1e944a94c496e302073435/server/bin/remote-cli:/u/morenodu/.pixi/bin:/u/morenodu/.pixi/bin:/opt/apps/miniconda/py312_24.7.1-0/bin:/opt/apps/miniconda/py312_24.7.1-0/condabin:/usr/share/Modules/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/u/morenodu/.vscode-server/extensions/ms-python.debugpy-2025.8.0-linux-x64/bundled/scripts/noConfigScripts:/u/morenodu/.vscode-server/data/User/globalStorage/github.copilot-chat/debugCommand"
export CONDA_PREFIX="/u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow"
export PIXI_PROJECT_MANIFEST="/u/morenodu/git_repos/compound-flooding-tropical-cyclones/pixi.toml"
export PIXI_IN_SHELL="1"
export PIXI_PROJECT_ROOT="/u/morenodu/git_repos/compound-flooding-tropical-cyclones"
export PIXI_PROJECT_VERSION="0.1.0"
export PIXI_EXE="/opt/apps/pixi/pixi"
export PIXI_PROJECT_NAME="COMPASS_UC3"
export CONDA_DEFAULT_ENV="COMPASS_UC3:compass-wflow"
export PIXI_ENVIRONMENT_NAME="compass-wflow"
export PIXI_ENVIRONMENT_PLATFORMS="linux-64,win-64"
export PIXI_PROMPT="(COMPASS_UC3:compass-wflow) "
. "/u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/activate.d/gdal-activate.sh"
. "/u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/activate.d/geotiff-activate.sh"
. "/u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/activate.d/libarrow_activate.sh"
. "/u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/activate.d/libglib_activate.sh"
. "/u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/activate.d/libxml2_activate.sh"
. "/u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/activate.d/proj4-activate.sh"
source /u/morenodu/git_repos/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/share/bash-completion/completions/*

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
