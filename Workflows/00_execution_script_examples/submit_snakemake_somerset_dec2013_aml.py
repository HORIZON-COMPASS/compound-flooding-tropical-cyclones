"""
Submit the Somerset Dec-2013 Snakemake workflow to the AzureML `eloise-cluster`.

This is the AzureML equivalent of
    00_execution_script_examples/run_snakemake_all_workflow_somerset_dec2013.sh

The heavy work runs on the `eloise-cluster` AmlCompute cluster; this script only submits.
It can therefore be run from anywhere that can reach the workspace.

From a compute instance:
    python Workflows/00_execution_script_examples/submit_snakemake_somerset_dec2013_aml.py

From a local machine:
    pip install azure-ai-ml azure-identity
    az login --tenant 17f18161-20d7-4746-87fd-50fe3e3b6619
    python Workflows/00_execution_script_examples/submit_snakemake_somerset_dec2013_aml.py

Workspace facts this script is wired to (queried from the workspace):
    subscription     d1d42766-9d6a-4571-8b38-bbc3e60604bb
    resource group   rg-ukcrcollab-ukcrcompassflood
    workspace        mlw-ukcrcompassflood-uksouth-01
    cluster          eloise-cluster  (Standard_E16ds_v4 -> 16 vCPU, 128 GB RAM)
    file datastore   workspaceworkingdirectory  (backs /home/azureuser/cloudfiles/code)
"""

from pathlib import Path

from azure.ai.ml import MLClient, Output, command
from azure.ai.ml.entities import BuildContext, Environment
from azure.identity import (
    ChainedTokenCredential,
    DefaultAzureCredential,
    InteractiveBrowserCredential,
    ManagedIdentityCredential,
)

# --------------------------------------------------------------------------------------
# 0. Paths
# --------------------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]       # .../hydromt_v1_somerset
CODE_DIR = REPO_ROOT / "Workflows"                    # snapshot uploaded with the job
# Build context: only the Dockerfile + conda yml. AzureML hashes this whole folder, so any
# other file in it (e.g. this script) would trigger a full image rebuild when edited.
DOCKER_DIR = Path(__file__).parent / "aml_image"

CONFIG_FILE = "config_general_somerset_TESTS.yml"
SNAKEFILE = "snakefile_all_wflow_sfincs.smk"
N_CORES = 16                                          # Standard_E16ds_v4

# --------------------------------------------------------------------------------------
# 1. Client
#    On a compute instance the managed identity works; on a laptop it does not exist, so
#    fall back to the az-CLI / VS Code / env-var chain and finally an interactive browser
#    login. Order matters: ManagedIdentityCredential is cheap to fail off-Azure.
# --------------------------------------------------------------------------------------
TENANT_ID = "17f18161-20d7-4746-87fd-50fe3e3b6619"

credential = ChainedTokenCredential(
    ManagedIdentityCredential(client_id="368335db-596d-437c-a2de-0c3dda300b73"),
    DefaultAzureCredential(exclude_managed_identity_credential=True),
    InteractiveBrowserCredential(tenant_id=TENANT_ID),
)

ml_client = MLClient(
    credential=credential,
    subscription_id="d1d42766-9d6a-4571-8b38-bbc3e60604bb",
    resource_group_name="rg-ukcrcollab-ukcrcompassflood",
    workspace_name="mlw-ukcrcompassflood-uksouth-01",
)

# --------------------------------------------------------------------------------------
# 2. Environment
#    Built from aml_image.Dockerfile rather than a bare conda spec: the workflow shells
#    out to Julia/Wflow.jl and to apptainer for the SFINCS .sif, neither of which conda
#    provides. AzureML builds the image on first submit and caches it; later submits with
#    an unchanged build context reuse it.
# --------------------------------------------------------------------------------------
compass_v1_env = Environment(
    name="compass-v1",
    description="HydroMT v1 + Wflow (julia) + SFINCS (apptainer) + Snakemake, pinned to pixi compass-v1.",
    build=BuildContext(path=str(DOCKER_DIR), dockerfile_path="aml_image.Dockerfile"),
)

# --------------------------------------------------------------------------------------
# 3. The command
#    `root_dir` in config_general_somerset_dec2013.yml points at the file share
#    (/home/azureuser/cloudfiles/code/Users/eloise.matthews). That path does not exist on
#    a cluster node, so we mount the same folder from the `workspaceworkingdirectory`
#    datastore and override root_dir on the Snakemake CLI (`--config root_dir=...`),
#    which the snakefiles honour (snakefile_all_wflow_sfincs.smk: root_dir = config['root_dir']).
#
#    The data catalogs (03_data_catalogs/data_catalog_MO_azure_UKCP.yml) hold ~18 absolute
#    `uri:` paths under that same share path, and hydromt resolves them verbatim. Rather
#    than rewrite the catalogs (the compute-instance runs depend on those paths), symlink
#    the expected location to the mount so both forms resolve to the same files.
# --------------------------------------------------------------------------------------
SHARE_PATH = "/home/azureuser/cloudfiles/code/Users/eloise.matthews"

snakemake_cmd = (
    "set -euo pipefail && "
    f"mkdir -p $(dirname {SHARE_PATH}) && "
    f"ln -sfn ${{{{outputs.workdir}}}} {SHARE_PATH} && "
    # snakefile_wflow.smk hardcodes --project=~/.julia/environments/wflow1, and julia
    # expands the tilde, so point that at the project baked into the image.
    'mkdir -p "$HOME/.julia/environments" && '
    'ln -sfn /opt/julia-depot/environments/wflow1 "$HOME/.julia/environments/wflow1" && '
    "cd 02_workflow_rules && "
    "echo 'root_dir -> ${{outputs.workdir}}' && "
    "snakemake --unlock "
    "  -s ${{inputs.snakefile}} "
    "  --configfile ../01_config_snakemake/${{inputs.config_file}} "
    "  --config root_dir=${{outputs.workdir}} && "
    "snakemake "
    "  -s ${{inputs.snakefile}} "
    "  --configfile ../01_config_snakemake/${{inputs.config_file}} "
    "  --config root_dir=${{outputs.workdir}} "
    "  --cores ${{inputs.n_cores}} "
    "  --latency-wait 120 "
    "  --wait-for-files "
    "  --forceall "
    "  --printshellcmds "
    "  --rerun-incomplete"
)

command_job = command(
    display_name="cluster_test",
    experiment_name="compass-somerset-tests",
    description="Snakemake all_wflow_sfincs workflow, Somerset Dec-2013 compound event (CF0).",
    code=str(CODE_DIR),
    command=snakemake_cmd,
    environment=compass_v1_env,
    inputs={
        "snakefile": SNAKEFILE,
        "config_file": CONFIG_FILE,
        "n_cores": N_CORES,
    },
    outputs={
        # Read-write mount of Users/eloise.matthews so the workflow sees 01_Data,
        # 02_Models (incl. 00_executables) and writes 03_Runs straight back to the share.
        "workdir": Output(
            type="uri_folder",
            path="azureml://datastores/workspaceworkingdirectory/paths/Users/eloise.matthews",
            mode="rw_mount",
        ),
    },
    environment_variables={
        "HDF5_USE_FILE_LOCKING": "FALSE",
        "PYTHONNOUSERSITE": "1",
        "JULIA_DEPOT_PATH": "/opt/julia-depot",
        # Keep numeric libs from oversubscribing the 16 cores Snakemake is scheduling on.
        "OMP_NUM_THREADS": "1",
        "NUMEXPR_MAX_THREADS": str(N_CORES),
    },
    compute="eloise-cluster",
    instance_count=1,
    shm_size="8g",
)

if __name__ == "__main__":
    returned_job = ml_client.jobs.create_or_update(command_job)
    print(f"Submitted: {returned_job.name}")
    print(f"Studio:    {returned_job.studio_url}")
    ml_client.jobs.stream(returned_job.name)
