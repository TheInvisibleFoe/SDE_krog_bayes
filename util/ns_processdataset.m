function [results]=ns_processdataset(obs,models,misc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine runs the nested sampling algorithm on a series of observed displacements
% (obs) for a range of models.
% For each model, the percentiles, means, and standard deviations are calculated
% via the 'ns_analyze' routine.
% Conclusively, the information yield in the data set is compared to that of replicated
% data for the same model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic
% Clean up data
obs = obs(:,~isnan(sum(obs,1)));


% Run nested sampling algorithm for each model
parfor i=1:length(models);
    [results(i).logZ,results(i).H,results(i).samples]...
        =ns_algorithm(obs,models(i));
end

% Calculate total evidence for the models
logZ_tot = log(0);
for i=1:length(models);
    logZ_tot = ns_logsumexp2(logZ_tot,results(i).logZ(1));
end

% Calculate normalized evidence for the models and more
for i = 1:length(models);
    results(i).Z_norm = exp(results(i).logZ(1) - logZ_tot);
    results(i).logZ_error = sqrt(results(i).H(1)/models(i).options.nwalkers);
    [results(i).percentiles,results(i).param_mean,results(i).param_stddev,results(i).maxLpar]...
        =ns_analyze(results(i).samples,models(i),misc);
end

% Calculate probability of  H(replicated obs) > H(obs)
for i = 1:length(models)
   results(i).prob = ns_infcheck(obs,models(i),results(i).logZ,results(i).samples);
end

%Print a summary of the results to a text file if wanted
if isfield(misc,'nssummary')
  ns_print(results,models,misc)
end

disp('Data processing complete');
toc

