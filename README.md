# Risk Loss aversion task and modelling
For any question please contact Caroline Charpentier (ccharpen@caltech.edu)

## Risk and Loss aversion task
The 'task' folder contains Matlab code for running the risk and loss aversion task (~10-15 minutes) using Cogent for stimuli presentation. Note that Cogent is not supported anymore, but can still be downloaded here: http://www.vislab.ucl.ac.uk/cogent_2000.php.
Therefore, I would recommend recoding this task anyway using PsychoPy (https://www.psychopy.org/) for in-lab studies or jsPsych (https://www.jspsych.org/) for online studies.

The task code ('Main_task_script.m') also contains a rough estimation of the risk and loss aversion parameters, although the code included in the 'modelling' folder is more up-to-date in terms of better optimization of the modelling and parameter estimation procedure.

Please refer to the 'Task_and_modelling_description.doc' in the 'task' folder for further instructions and details.

## Risk and Loss aversion modelling
The 'modelling' folder contains the Matlab code necessary to estimate risk and loss aversion from a multi-session gambling task, used for modelling the data of Stamatis, Puccetti et al (2020) Behaviour Research and Therapy (https://www.sciencedirect.com/science/article/abs/pii/S0005796720300607?via%3Dihub)

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
