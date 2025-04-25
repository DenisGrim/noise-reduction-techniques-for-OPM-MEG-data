
folder of vbmeg and raw data might has to be adjusted (former in make\*.m and latter in set_parameters.m)

the 'figures' folder contains a copy of each figure-folder created by the pipeline for each noise correction technique
To remake it, four runs of make_all are needed (fourth technique is 'none'), each run with a different technique chosen in set_parameters.m.


most of the code is an exact copy of https://vbmeg.atr.jp/docs/v30/static/vbmeg3_opm_processing.html

edited/added files are (all in /program):
make_all.m
make_one_subject_one_task.m
set_parameters.m
den_spm_opm_vslm.m (exact copy of https://github.com/spm/spm/blob/22cab57139a07650870bfd09cc76077a7688c63d/toolbox/MEEGtools/spm_opm_vslm.m#L4 other than name)
toolbox/get_SNR.m
visualization/show_trial_average.m
preprocess/apply_car.m
preprocess/appy_harmonicfc.m
