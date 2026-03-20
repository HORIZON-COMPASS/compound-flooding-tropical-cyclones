. "/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/deactivate.d/gdal-deactivate.sh"
. "/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/deactivate.d/geotiff-deactivate.sh"
. "/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/deactivate.d/libglib_deactivate.sh"
. "/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/deactivate.d/libxml2_deactivate.sh"
. "/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/deactivate.d/proj4-deactivate.sh"
. "/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/deactivate.d/udunits2-deactivate.sh"
export PATH="/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/bin:/home/users/gregory.munday/.pixi/bin:/home/users/gregory.munday/.pixi/bin:/home/users/gregory.munday/.juliaup/bin:/home/users/gregory.munday/.pixi/bin:/home/users/gregory.munday/.local/bin:/home/users/gregory.munday/bin:/opt/moose-client-wrapper/bin:/opt/conda/condabin:/usr/local/bin:/usr/bin:/bin:/sbin:/usr/local/sbin:/usr/sbin:/opt/ukmo/idl/ukmo/bin:/data/apps/maple/2021/bin:/opt/cycle/jetpack/bin"
export CONDA_PREFIX="/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow"
export PIXI_PROJECT_NAME="COMPASS_UC3"
export PIXI_PROJECT_MANIFEST="/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/pixi.toml"
export PIXI_PROJECT_VERSION="0.1.0"
export PIXI_IN_SHELL="1"
export PIXI_PROJECT_ROOT="/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones"
export PIXI_EXE="/home/users/gregory.munday/.pixi/bin/pixi"
export CONDA_DEFAULT_ENV="COMPASS_UC3:compass-wflow"
export PIXI_ENVIRONMENT_NAME="compass-wflow"
export PIXI_ENVIRONMENT_PLATFORMS="linux-64,win-64"
export PIXI_PROMPT="(COMPASS_UC3:compass-wflow) "
. "/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/activate.d/gdal-activate.sh"
. "/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/activate.d/geotiff-activate.sh"
. "/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/activate.d/libarrow_activate.sh"
. "/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/activate.d/libglib_activate.sh"
. "/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/activate.d/libxml2_activate.sh"
. "/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/activate.d/proj4-activate.sh"
. "/home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/etc/conda/activate.d/udunits2-activate.sh"
source /home/users/gregory.munday/Documents/compound-flooding-tropical-cyclones/.pixi/envs/compass-wflow/share/bash-completion/completions/*

# shellcheck shell=bash
pixi() {
    local first_arg="${1-}"

    "${PIXI_EXE-}" "$@" || return $?

    case "${first_arg-}" in
        add|a|remove|rm|install|i)
            eval "$("$PIXI_EXE" shell-hook --change-ps1 false)"
            hash -r
            ;;
    esac || :

    return 0
}

export PS1="(COMPASS_UC3:compass-wflow) ${PS1:-}"
