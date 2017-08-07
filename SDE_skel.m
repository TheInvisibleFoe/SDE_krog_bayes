function SDE_skel()

% This skeleton file loads the data specified by the 'misc.data_id' parameter
% and initiates Bayesian data analysis via the Nested sampling algorithm.
%
% The prior ranges stored in the 'ranges' array specify the user set limit on the
% model parameters to be sampled from.
%
% The 'option' struct holds various parameters controlling the Nested sampling routine 
% 
% The 'model' struct specify the model dependent functions used for calculating the likelihood
% and converting [0;1] numbers to model parameters as well as data scaling. 
%
% NB:
% For the 'SDE' library, these function are stored in the 'util' directory so the command
%      "addpath('util')"
% must be called in MatLab in order to allow running of the 'SDE_skel.m' routine.
 

%Specify location of data
misc.data_id = 'samples/ftrack';

data_path = [misc.data_id,'.txt'];

%load data 

data = importdata(data_path);

%Specify prior ranges
Dmin=10^(-4);
Dmax=10^2;   

amin = -2;
amax = 2;

fmin = -1;
fmax = 1;

noisemin=0; 
noisemax=100;

% ranges - a 2x4 array of minimum and maximum values for the 4 parameters:
%   - Mobility coefficient, mobility exponent, force, and measurement error
ranges=[Dmin Dmax ; amin amax; fmin fmax; noisemin noisemax];

Tmin = 10; % Required minimum length of a scaled data set
nmax = floor((length(data)-1) / Tmin);
if ~(nmax > 1)
   fprintf('Data set too short to be scaled to length >= %i',Tmin);
elseif nmax > 10
   nmax = 10;
end


%Specify options
options.nwalkers=300;   % Number of walkers to be generated
options.stoprat=10^(-6);% Ratio to stop sampling
options.nsteps=30;      % Number of walkers in parameter space
options.nmax = nmax;    % Maximum time scaling factor for data
options.trackmax = 100; % Number of replicated tracks to compare information for

final=@(x) x(end); % extract last element in a vector

%Specify the models
MM = [0 0; 1 0; 0 1; 1 1];
n_perm = 2; % Number of permanently active parameters
for i=1:size(MM,1);
  options.lengthu=sum(MM(i,:))+n_perm;
  models(i).options=options;
  models(i).logl=@(obs,theta) SDE_logl(obs,SDE_params(theta,MM(i,:)));
  models(i).logl_n=@(obs,theta,n) SDE_logl_m(obs,SDE_params(theta,MM(i,:)),n);
  models(i).scaling =@(obs,n) SDE_scaling(obs,n); % Function that scales data
  models(i).replicate =@(obs,theta,n) SDE_replicate(obs,SDE_params(theta,MM(i,:)),n);
  models(i).invprior=@(u) SDE_invprior(u,ranges,MM(i,:));
  models(i).labels=[1:n_perm];
  for j=1:length(MM(i,:));
    if MM(i,j)==1
      models(i).labels=[models(i).labels j+n_perm];
    end
  end
end

%Percentiles to be calculated
misc.percentiles_at=[0.02 0.16 0.5 0.84 0.98];

%Labels for the parameters
misc.labels=...
['Step deviation: ';...
 'Exponent:       ';...
 'Force:          ';...
 'Measurem. err.: '];

misc.titles=...
{'Percentile at:  ',...
 ' MaxL@',...
 ' Mean	+/-  dev.'};

%Tell ns_print to write a summary-file
misc.nssummary=['_results.txt'];

%Run the nested sampling algorithm for all models and compile results
[results] = ns_processdataset(data,models,misc);

path=[misc.data_id,'_output'];

%Save results
save(path,'results')
save(path,'data','-append')


