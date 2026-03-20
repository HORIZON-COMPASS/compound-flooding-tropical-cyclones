rule add_forcing_coastal_meteo_sfincs:
    input:
        msk_file = join(root_dir, dir_models, "{region}", "{runname}", "sfincs", "sfincs.msk"),
    params:
        dir_run_no_forcing = directory(join(root_dir, dir_models, "{region}", "{runname}", "sfincs")),
        dir_run_with_forcing = directory(join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}")),
        data_cats = get_datacatalog,
        wind_forcing = get_wind_forcing,
        start_time = get_starttime, 
        end_time = get_endtime,
        use_dfm = get_use_dfm,
        coastal_ts = get_coastal_ts,
        dfm_output = lambda wildcards: "dfm_output_event_"+ config['runname_ids'][wildcards.runname]["dfm_res"] + "_" + config['runname_ids'][wildcards.runname]["bathy"] + "_" + config['runname_ids'][wildcards.runname]["tidemodel"] + "_" + config['runname_ids'][wildcards.runname]["wind_forcing"],
        utmzone = get_utmzone,
    output:
        bzs_file = join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}", "sfincs.bzs"),
    script:
        join( '..', "04_scripts", "model_building", "sfincs", "update_sfincs_coastal_forcing.py")

rule update_dis_forcing_sfincs:
    input:
        bzs_file = join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}", "sfincs.bzs"),
        wflow_output = join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}", "events", "run_default", "output_scalar.nc"),
    params:
        dir_run_with_forcing = directory(join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}")),
        wflow_root_forcing = directory(join(root_dir, dir_runs, "{region}", "{runname}", "wflow","event_precip_{forcing}")),
        data_cats = get_datacatalog,
    output:
        dis_file = join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}", "sfincs.dis"),
    script:
        join( '..', "04_scripts", "model_building", "sfincs", "update_sfincs_dis_forcing.py")


rule run_sfincs_model:
    input:
#        batchfile = "{dir_run}"+"/sfincs_"+"{runname}"+"/run_sfincs.bat"
        dis_file = join(root_dir,  dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}", "sfincs.dis"),
    params:
        dir_run_with_forcing = directory(join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}")),
        exe = join(root_dir, dir_models, "00_executables", "SFINCS_v2.1.1_Dollerup_release_exe", 'sfincs.exe'),
        currentdir = curdir
    output:
        mapout = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}", "sfincs_map.nc"),
    run:
        if os.name == 'nt':
            import subprocess
            print(f"Running Sfincs model at {params.dir_run_with_forcing} with bin {params.exe}")
            print("Executing SFINCS...")
            with open (join(params.dir_run_with_forcing,"sfincs.log"), "w") as f:
                subprocess.run([str(params.exe)], stdout=f, cwd=params.dir_run_with_forcing)
                print("Finished running")
        if os.name == 'posix': 
            # Met Office uses Apptainer instead of Docker 
            # shell("apptainer pull docker://deltares/sfincs-cpu")
            shell("apptainer run {params.dir_run_with_forcing}/sfincs-cpu_latest.sif")


rule sfincs_plot_floodmap:
    input:
        mapout = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}", "sfincs_map.nc"),
    params:
        dir_run = directory(join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}")),
        dir_model_no_forcing = directory(join(root_dir, dir_models, "{region}", "{runname}", "sfincs")),
        datacat = get_datacatalog
    output:
        figure = join(root_dir, dir_runs, "{region}", "{runname}", "sfincs","event_precip_{forcing}", "plot_output", "sfincs_basemap.png")  
    script:
        join(curdir,  '..', "04_scripts", "postprocessing", "sfincs", "sfincs_postprocess.py")
