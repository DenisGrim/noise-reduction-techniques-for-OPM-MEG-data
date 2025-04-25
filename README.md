
folder of vbmeg and raw data might has to be adjusted (former in make(...).m and latter in set_parameters.m)

the 'figures' folder contains a copy of each figure-folder created by the pipeline for each noise correction technique
To remake it, four runs of make_all are needed (fourth technique is 'none'), each run with a different technique chosen in set_parameters.m.


most of the code is an exact copy of https://vbmeg.atr.jp/docs/v30/static/vbmeg3_opm_processing.html

