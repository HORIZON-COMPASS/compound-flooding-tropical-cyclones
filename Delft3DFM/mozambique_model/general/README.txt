This folder is for:
 - general files needed to run the model which do not change in between runs.


base_model_settings.mdu:
 - Based on GTSM model settings, removed some keywords that caused errors with Hydrolib, most importantly the "SelfAttractionLoading" (was set to 1). For a local model this is likely not a big issue.