%This script loads individual data files from a risk/loss aversion task
%that was run 3 times for each participant: at baseline, during
%anticipation of a stressor, and during recovery (after the stressor).
%It then collapse the data of all three sessions together and estimates two
%models:
%- a model where loss aversion (lambda), risk preference (rho) and choice
%consistency (mu) are estimated across all trials (calling functions
%LL_function_alldata.m and generate_choice_alldata.m)
%- a model where each of the 3 parameters are estimated separately for each
%session (baseline, anticipation, recovery), leading to 9 parameters total.
%This model calls LL_function_split_per_session.m and
%generate_choice_split_per_session.m
%In order for the gambling data to be passed into the modelling function,
%it needs to have the following format: matrix "P", each row is a trial
% column 2 is the value of the sure option, 
% column 3 is the value of the potential loss,
% column 4 is the value of the potential gain, 
% column 8 is the binary choice variable: 1 if gamble, 0 if sure option.
% column 10 here depicts the session (baseline=1, anticipation=2,
% recovery=3), necessary for Model 2

%------ Caroline Charpentier (ccharpen@caltech.edu) November 2018 --------

clear all
close all

fs = filesep;

data_dir = 'C:\Users\Caroline\Box Sync\Gambling_Stress_data\Data'; %change with directory where data (subject folders) are saved
cd(data_dir)
folders = dir([data_dir fs 'Sub*']);
sub_list = {folders.name};

% init the randomization screen
RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));

%Setup solver opt 
% for fminunc (Jeff's settings)
npar = 5; %number of parameters
opts = optimoptions(@fminunc, ...
        'Algorithm', 'quasi-newton', ...
        'Display','off', ...
        'MaxFunEvals', 50 * npar, ...
        'MaxIter', 50 * npar,...
        'TolFun', 0.01, ...
        'TolX', 0.01);
    
sim_nb = 10; %number of simulations for recovery (for reliable results, this should be increased to 50 or 100)

fitResult = struct();

Model1_Params          = zeros(length(sub_list),3);
Model1_Loglikelihood   = zeros(length(sub_list),1);
Model1_ExitFlag        = zeros(length(sub_list),1);
Model1_BIC             = zeros(length(sub_list),1);
Model1_AIC             = zeros(length(sub_list),1);
Model1_PseudoR2        = zeros(length(sub_list),1);
Model1_RecoveredParams = zeros(length(sub_list),3);
Model1_RecoveredStd    = zeros(length(sub_list),3);

Model2_Params          = zeros(length(sub_list),9);
Model2_Loglikelihood   = zeros(length(sub_list),1);
Model2_ExitFlag        = zeros(length(sub_list),1);
Model2_BIC             = zeros(length(sub_list),1);
Model2_AIC             = zeros(length(sub_list),1);
Model2_PseudoR2        = zeros(length(sub_list),1);
Model2_RecoveredParams = zeros(length(sub_list),9);
Model2_RecoveredStd    = zeros(length(sub_list),9);

parfor i=1:length(sub_list)
    
    disp(sub_list{i})
    sub_dir = [data_dir fs sub_list{i}];
    
    %create a big matrix with trial data from the 3 timepoints (baseline,
    %anticipation and recovery)
    fname = dir([sub_dir fs 'base_*']);
    DATA = load([sub_dir fs fname.name],'DATA');
    data_base = DATA.DATA.TrialsCombined2;
    
    fname = dir([sub_dir fs 'anti_*']);
    DATA = load([sub_dir fs fname.name],'DATA');
    data_anti = DATA.DATA.TrialsCombined2;
    
    fname = dir([sub_dir fs 'recov_*']);
    DATA = load([sub_dir fs fname.name],'DATA');
    data_recov = DATA.DATA.TrialsCombined2;
    
    data_all = [data_base ones(length(data_base(:,1)),1); data_anti 2*ones(length(data_anti(:,1)),1); ...
        data_recov 3*ones(length(data_recov(:,1)),1)];
    P = data_all; %for models
    
    tr_nb = length(P(:,1));
    
    n_good_trials(i,:) = [length(data_base(:,1)) length(data_anti(:,1)) length(data_recov(:,1)) tr_nb];
    
    %% Fit Model 1
    %global model: 1 mu, 1 rho, 1 lambda across all data (excluding staircase)
    disp([sub_list{i} ' - Fitting Mod 1'])
    npar = 3;
    params_rand=[]; 
    for i_rand=1:20*npar
        start = randn(1,npar); %also possible: sample from a gamma distribution (rather than normal)
        [paramtracker, lltracker, ex] = fminunc(@LL_function_alldata, start, opts, P);
        params_rand=[params_rand; paramtracker lltracker ex];
    end
    bad_fit = params_rand(:,npar+2)<=0;
    if sum(bad_fit)<20*npar
        params_rand(bad_fit,:)=[];
    end
    [~,ids] = sort(params_rand(:,npar+1));
    best_params = params_rand(ids(1),:);
    best_params(1) = exp(best_params(1));   % mu [0 Inf]
    best_params(2) = exp(best_params(2));   % lambda [0 Inf]
    best_params(3) = exp(best_params(3));   % rho [0 Inf]
    BIC = 2*best_params(npar+1) + npar*log(tr_nb);
    AIC = 2*best_params(npar+1) + 2*npar;
    pseudoR2 = 1 + best_params(npar+1)/(log(0.5)*tr_nb);
    
    Model1_Params(i,:)        = best_params(1:npar);
    Model1_Loglikelihood(i,1) = best_params(npar+1);
    Model1_ExitFlag(i,1)      = best_params(npar+2);
    Model1_BIC(i,1)           = BIC;
    Model1_AIC(i,1)           = AIC;
    Model1_PseudoR2(i,1)      = pseudoR2;
    
    disp([sub_list{i} ' - Recov Mod 1'])
    params = best_params(1:npar);
    param_r = zeros(sim_nb,npar); %initiate list of recovered parameter
    for s=1:sim_nb
        P_pred = generate_choice_alldata(params,P); %generate choice set
        gen_ch = P_pred(:,4); %set of generated choices, column nb may need to be changed!
        P_new = P;
        P_new(:,8) = gen_ch;
        %fit model to generated choice set
        params_rand=[]; 
        for i_rand=1:20*npar
            start = randn(1,npar); %also possible: sample from a gamma distribution (rather than normal)
            [paramtracker, lltracker, ex] = fminunc(@LL_function_alldata, start, opts, P_new);
            params_rand=[params_rand; paramtracker lltracker ex];
        end
        bad_fit = params_rand(:,npar+2)<=0;
        if sum(bad_fit)<20*npar
            params_rand(bad_fit,:)=[];
        end
        [~,ids] = sort(params_rand(:,npar+1));
        best_params = params_rand(ids(1),:);
        best_params(1) = exp(best_params(1));   % mu [0 Inf]
        best_params(2) = exp(best_params(2));   % lambda [0 Inf]
        best_params(3) = exp(best_params(3));   % rho [0 Inf]
        param_r(s,:) = best_params(1:npar);
    end
    Model1_RecoveredParams(i,:) = mean(param_r,1);
    Model1_RecoveredStd(i,:) = std(param_r,1);
    
    %% Fit Model 2
    %model where each parameter is estimated separately for each session
    %(baseline, anticipation & recovery)
    disp([sub_list{i} ' - Fitting Mod 2'])
    npar = 9;
    params_rand=[]; 
    for i_rand=1:20*npar
        start = randn(1,npar); %also possible: sample from a gamma distribution (rather than normal)
        [paramtracker, lltracker, ex] = fminunc(@LL_function_split_per_session, start, opts, P);
        params_rand=[params_rand; paramtracker lltracker ex];
    end
    bad_fit = params_rand(:,npar+2)<=0;
    if sum(bad_fit)<20*npar
        params_rand(bad_fit,:)=[];
    end
    [~,ids] = sort(params_rand(:,npar+1));
    best_params = params_rand(ids(1),:);
    for p = 1:npar
        best_params(p) = exp(best_params(p));   % transform to [0 Inf]
    end
    BIC = 2*best_params(npar+1) + npar*log(tr_nb);
    AIC = 2*best_params(npar+1) + 2*npar;
    pseudoR2 = 1 + best_params(npar+1)/(log(0.5)*tr_nb);
    
    Model2_Params(i,:)        = best_params(1:npar);
    Model2_Loglikelihood(i,1) = best_params(npar+1);
    Model2_ExitFlag(i,1)      = best_params(npar+2);
    Model2_BIC(i,1)           = BIC;
    Model2_AIC(i,1)           = AIC;
    Model2_PseudoR2(i,1)      = pseudoR2;
    
    disp([sub_list{i} ' - Recov Mod 2'])
    params = best_params(1:npar);
    param_r = zeros(sim_nb,npar); %initiate list of recovered parameter
    for s=1:sim_nb
        P_pred = generate_choice_split_per_session(params,P); %generate choice set
        gen_ch = P_pred(:,4); %set of generated choices, column nb may need to be changed!
        P_new = P;
        P_new(:,8) = gen_ch;
        %fit model to generated choice set
        params_rand=[]; 
        for i_rand=1:20*npar
            start = randn(1,npar); %also possible: sample from a gamma distribution (rather than normal)
            [paramtracker, lltracker, ex] = fminunc(@LL_function_split_per_session, start, opts, P_new);
            params_rand=[params_rand; paramtracker lltracker ex];
        end
        bad_fit = params_rand(:,npar+2)<=0;
        if sum(bad_fit)<20*npar
            params_rand(bad_fit,:)=[];
        end
        [~,ids] = sort(params_rand(:,npar+1));
        best_params = params_rand(ids(1),:);
        for p = 1:npar
            best_params(p) = exp(best_params(p));   % transform to [0 Inf]
        end
        param_r(s,:) = best_params(1:npar);
    end
    Model2_RecoveredParams(i,:) = mean(param_r,1);
    Model2_RecoveredStd(i,:) = std(param_r,1);
    
end

fitResult.Model1.Params          = Model1_Params;
fitResult.Model1.Loglikelihood   = Model1_Loglikelihood;
fitResult.Model1.ExitFlag        = Model1_ExitFlag;
fitResult.Model1.BIC             = Model1_BIC;
fitResult.Model1.AIC             = Model1_AIC;
fitResult.Model1.PseudoR2        = Model1_PseudoR2;
fitResult.Model1.Recovery.MeanParams   = Model1_RecoveredParams;
fitResult.Model1.Recovery.StdParams    = Model1_RecoveredStd;

fitResult.Model2.Params          = Model2_Params;
fitResult.Model2.Loglikelihood   = Model2_Loglikelihood;
fitResult.Model2.ExitFlag        = Model2_ExitFlag;
fitResult.Model2.BIC             = Model2_BIC;
fitResult.Model2.AIC             = Model2_AIC;
fitResult.Model2.PseudoR2        = Model2_PseudoR2;
fitResult.Model2.Recovery.MeanParams   = Model2_RecoveredParams;
fitResult.Model2.Recovery.StdParams    = Model2_RecoveredStd;

save('Results_models_nobounds_Nov2018.mat','fitResult')