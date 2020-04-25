# Risk Loss aversion task and modelling

### Risk and Loss aversion task

### Risk and Loss aversion modelling
The modelling folder contains the Matlab code necessary to estimate risk and loss aversion from a multi-session gambling task.

The modelling script to be run is the "run_models_nobounds_Nov2018.m", which calls the other functions.

This script first loads individual data files from a risk/loss aversion task that was run 3 times for each participant: at baseline, during anticipation of a stressor, and during recovery (after the stressor). See files in SubExample folder for an example of raw data files for one subject.

It then collapse the data of all three sessions together and estimates two models:
- a model where loss aversion (lambda), risk preference (rho) and choice consistency (mu) are estimated across all trials (Model 1, calling functions LL_function_alldata.m and generate_choice_alldata.m)
- a model where each of the 3 parameters are estimated separately for each session (baseline, anticipation, recovery), leading to 9 parameters total. This model (Model 2) calls LL_function_split_per_session.m and generate_choice_split_per_session.m

In order for the gambling data to be passed into the modelling function, it needs to have the following format: matrix "P", each row is a trial and:
- column 2 is the value of the sure option,
- column 3 is the value of the potential loss,
- column 4 is the value of the potential gain,
- column 8 is the binary choice variable: 1 if gamble, 0 if sure option.
- column 10 here depicts the session (baseline=1, anticipation=2, recovery=3), necessary for Model 2

The script's output file is saved as "Results_models_nobounds_Sept2018.mat".

For any question please contact Caroline Charpentier (ccharpen@caltech.edu)
