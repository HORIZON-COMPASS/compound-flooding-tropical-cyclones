import argparse
from pathlib import Path
import yaml


def update_config(experiment_name: str, forcing_dataset: str, event_id: str) -> None:
    with open("Workflows/01_config_snakemake/config_general_somerset.yml", "r") as file:
        yml = yaml.safe_load(file)

    # Update the config
    original_key = list(dict(yml['runname_ids']).keys())[0]
    exp_settings = dict(yml['runname_ids'][original_key]) # new dict structure
    exp_settings['forcing'] = f"{forcing_dataset}_{event_id}" # overwrite the forcing
    yml['runname_ids'][experiment_name[:4] + str(event_id)] = exp_settings  # add as new runname_id
    yml['runname_ids'].pop(original_key, None)

    with open(f"Workflows/01_config_snakemake/config_general_{forcing_dataset}_{event_id}.yml", "w") as file:
        yaml.dump(yml, file)

    output_dir = Path(yml['root_dir']) / yml['root_dir'] / yml['dir_runs'] / exp_settings['region'] / experiment_name / "sfincs" / f"event_precip_{forcing_dataset}_{event_id}"
    print(str(output_dir))


def update_data_catalog(forcing_dataset: str, event_id: str):
    with open("Workflows/03_data_catalogs/data_catalog_MO.yml", "r") as file:
        yml = yaml.safe_load(file)

    # Update the event_id in the relevant forcing entry
    template = dict(yml[forcing_dataset])
    path_base = yml[forcing_dataset]['meta']['source_url']
    path = path_base.replace('_id_', f'_{event_id}_')
    template['path'] = path
    yml[f"{forcing_dataset}_{event_id}"] = template

    with open("Workflows/03_data_catalogs/data_catalog_MO.yml", "w") as file:
        yaml.dump(yml, file)


def main(experiment_name: str, forcing_dataset: str, event_id: int) -> None:
    update_config(experiment_name, forcing_dataset, event_id)
    update_data_catalog(forcing_dataset, event_id)


if __name__ == "__main__":
    # Parse arguments from the batch script
    parser = argparse.ArgumentParser("update_config")
    parser.add_argument("--experiment", help="The name of the experiment you want to run", type=str)
    parser.add_argument("--dataset", help="The name of the forcing dataset you want to run with", type=str)
    parser.add_argument("--event", help="The event number related to the forcing dataset", type=int)
    args = parser.parse_args()
    main(args.experiment, args.dataset, args.event)