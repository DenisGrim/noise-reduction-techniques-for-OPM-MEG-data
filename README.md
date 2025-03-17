
folder of vbmeg and raw data might has to be adjusted (former in make(...).m and latter in set_parameters.m)

the 'figures' folder contains a copy of each figure-folder created by the pipeline for each noise correction technique
To remake it, three runs of make_all are needed, each run with a different function used in the noise correction step:
- apply_car -> common average referencing
- apply_harmonicfc -> harmonic field correction
- apply_hfc -> homogeneous field correction

and the figure-folder must be copied into a different folder after each run! Otherwise they will be overwritten


most of the code is an exact copy of https://vbmeg.atr.jp/docs/v30/static/vbmeg3_opm_processing.html

