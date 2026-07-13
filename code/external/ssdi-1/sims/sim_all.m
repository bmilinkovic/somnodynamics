%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up SS or VAR model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc




CON = tnet9x
varmod = exist('CON','var'); % for a VAR model, specify a connectivity matrix or a scalar dimension
if varmod % VAR model
	if isscalar(CON), n = CON; CON = ones(n); else, n = size(CON,1); end;
	defvar('r', 7   ); % VAR model order
	defvar('w', 1   ); % VAR coefficients decay parameter
else      % fully-connected state-space model
	defvar('n', 9   ); % microscopic dimension
	defvar('r', 3*n ); % hidden state dimension
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('moddir',   'results/for_thesis/'    ); % simulation directory
defvar('modname', 'sim_model' ); % model filename root
defvar('poptname', 'preopt_dd' );  % pre-optimisation filename root
defvar('optname',  'opt_dd'    ); % optimisation filename root
defvar('goptpname',  'goptp_dd'    ); % sub-optima distances filename root
defvar('goptoname',  'gopto_dd'    ); % optima distances filename root

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('rho',      0.9   ); % spectral norm (< 1)
defvar('rmii',     1     ); % residuals multiinformation; 0 for zero correlation
defvar('fres',     []    ); % frequency resolution (empty for automatic)
defvar('nsics',    0     ); % number of samples for spectral integration check (0 for no check)
defvar('mseed',    0     ); % model random seed (0 to use current rng state)
defvar('iseed',     0    ); % initialisation random seed (0 to use current rng state)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('nrunsp',    100        ); % pre-optimisation runs (restarts)
defvar('nitersp',   10000      ); % pre-optimisation iterations
defvar('gdesp',     2          ); % gradient-descent ES version (1 or 2)
defvar('gdsig0p',   1          ); % pre-optimisation (gradient descent) initial step size
defvar('gdlsp',     2          ); % gradient-descent "line search" parameters
defvar('gdtolp',    1e-10      ); % gradient descent convergence tolerance
defvar('histp',     true       ); % calculate optimisation history?
defvar('ppp',       false      ); % parallelise multiple runs?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('ctol',      1e-6       ); % hyperplane clustering tolerance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('niterso',   10000      ); % optimisation iterations
defvar('gdeso',     2          ); % gradient-descent ES version (1 or 2)
defvar('gdsig0o',   0.1        ); % optimisation (gradient descent) initial step size
defvar('gdlso',     2          ); % gradient-descent "line search" parameters
defvar('gdtolo',    1e-10      ); % gradient descent convergence tolerance
defvar('histo',     true       ); % calculate optimisation history?
defvar('ppo',       false      ); % parallelise multiple runs?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOR PLOTTING

defvar('gpterm',   'x11'       ); % Gnuplot terminal
defvar('gpscale',  [Inf,1.1]   ); % Gnuplot scale
defvar('gpfsize',  14          ); % Gnuplot font size
defvar('gpplot',   2           ); % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)




%%
% Generate random VAR or ISS model
%
% NOTE: parameters, etc., with 0 suffix are for untransformed model!

rstate = rng_seed(mseed);
V0 = corr_rand(n,rmii); % residuals covariance (actually correlation) matrix
if varmod
	mdescript = sprintf('%d-variable VAR(%d)',n,r);
	ARA0 = var_rand(CON,r,rho,w);           % random VAR model
	gc = var_to_pwcgc(ARA0,V0);             % causal graph
	[A0,C0,K0] = var_to_ss(ARA0);           % equivalent ISS model
	if isempty(fres)
		[fres,ierr] = var2fres(ARA0,V0);
	end
	modcomp = r*n*n + (n*(n+1))/2;          % model complexity
else
	mdescript = sprintf('%d-variable ISS(%d)',n,r);
	[A0,C0,K0] = iss_rand(n,r,rho);         % random ISS model
	gc = ss_to_pwcgc(A0,C0,K0,V0);          % causal graph
	if isempty(fres)
		[fres,ierr] = ss2fres(A0,C0,K0,V0);
	end
	modcomp = (2*n+1)*r + (n*(n+1))/2;      % model complexity
end
rng_restore(rstate);

% Spectral resolution accuracy check

% if nsics > 0
% 	assert(exist('mdim','var'),'For spectral accuracy check, must supply macro dimension ''mdim''');
% 	derr = dds_check(A,C,K,H,mdim,nsics); % spectral integration check
% 	if derr > 1e-12, fprintf(2,'WARNING: spectral DD calculation may be inaccurate!\n\n'); end
% end

% Model info

fprintf('\n---------------------------------------\n');
fprintf('Model            : %s\n',mdescript);
fprintf('---------------------------------------\n');
fprintf('Dimension        : %d\n',n);
fprintf('Complexity       : %d\n',modcomp);
fprintf('Freq. resolution : %d\n',fres);
fprintf('---------------------------------------\n\n');

% Get edgeweights for causal graph

eweight = gc/nanmax(gc(:)); 

% Save model

modfile = [fullfile(moddir,modname) '.mat'];
fprintf('*** saving model in ''%s''... ',modfile);
if varmod
	save(modfile,'mdescript','varmod','n','r','rho','rmii','V0','ARA0','A0','C0','K0','gc','fres','modcomp', 'eweight');
else
	save(modfile,'mdescript','varmod','n','r','rho','rmii','V0',       'A0','C0','K0','gc','fres','modcomp', 'eweight');
end
fprintf('done\n\n');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical dependence pre-optimisation (via proxy DD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Supply model for your data (see, e.g., sim_model.m).
%
% Specify a macroscopic dimension mdim and optimisation parameters, or accept
% defaults (see below). After running this script, run optimise_dd.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defvar('mdim',    8); % macroscopic dimension
assert(exist('mdim','var'),'Must supply macro dimension ''mdim''');

% Load model

modfile = [fullfile(moddir,modname) '.mat'];
fprintf('\n*** loading model from ''%s''... ',modfile);
load(modfile);
fprintf('done\n\n');

% Transform model to decorrelated and normalised form, and calculate CAK sequence

if varmod
	[ARA,V] = transform_var(ARA0,V0);       % transform model to decorrelated-residuals form
	[A,C,K] = var_to_ss(ARA);               % equivalent ISS model
	CAK = ARA;
else
	[A,C,K,V] = transform_ss(A0,C0,K0,V0);  % transform model to decorrelated-residuals form
	CAK = iss2cak(A,C,K);
end



fprintf('%s: pre-optimisation for m = %d\n\n',mdescript,mdim);

% Initialise optimisations

rstate = rng_seed(iseed);
L0p = rand_orthonormal(n,mdim,nrunsp); % initial (orthonormalised) random linear projections
rng_restore(rstate);

% Multiple optimisation runs

st = tic;
[doptp,Lp,convp,ioptp,soptp,cputp,ohistp] = opt_gd_ddx_mruns(CAK,L0p,nitersp,gdesp,gdsig0p,gdlsp,gdtolp,histp,ppp);
et = toc(st);

% Inverse-transform Lp back for un-decorrelated residuals

Loptp = itransform_subspace(Lp,V0);

fprintf('\noptimal dynamical dependence =\n'); disp(doptp');
fprintf('Simulation time = %s\n\n',datestr(seconds(et),'HH:MM:SS.FFF'));
fprintf('CPU secs per run = %7.4f +- %6.4f\n\n',mean(cputp),std(cputp));

% Get inter-optima subspace distances

goptp = gmetrics(Loptp);
goptpfile = fullfile(moddir,[goptpname '_mdim_' num2str(mdim) '_H.mat']);
save(goptpfile);

% Save pre-optimisation results

if histp
	poptfile = fullfile(moddir,[poptname '_mdim_' num2str(mdim) '_H.mat']);
else
	poptfile = fullfile(moddir,[poptname '_mdim_' num2str(mdim) '_N.mat']);
end
fprintf('\n*** saving pre-optimisation results in ''%s''... ',poptfile);
save(poptfile);
fprintf('done\n\n');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical dependence optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% You must run preoptimise_dd first!
%
% Specify a macroscopic dimension mdim and optimisation parameters, or accept
% defaults (see below).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate transfer function

if varmod
	H = var2trfun(ARA,fres);
else
	H = ss2trfun(A,C,K,fres);
end

fprintf('%s: optimisation (fres = %d) for m = %d\n',mdescript,fres,mdim);

% Cluster preoptimised hyperplanes

[uidx,~,nrunso] = Lcluster(goptp,ctol,doptp);

% Initialise optimisations

L0o = Lp(:,:,uidx);

% Multiple optimisation runs

st = tic;
[dopto,Lo,convp,iopto,sopto,cputo,ohisto] = opt_gd_dds_mruns(H,L0o,niterso,gdeso,gdsig0o,gdlso,gdtolo,histo,ppo);
et = toc(st);

% Inverse-transform Lo back for un-decorrelated residuals

Lopto = itransform_subspace(Lo,V0);

fprintf('\noptimal dynamical dependence =\n'); disp(dopto');
fprintf('Simulation time = %s\n\n',datestr(seconds(et),'HH:MM:SS.FFF'));
fprintf('CPU secs per run = %7.4f +- %6.4f\n\n',mean(cputo),std(cputo));

% Get inter-optima subspace distances

gopto = gmetrics(Lopto);
goptofile = fullfile(moddir,[goptoname '_mdim_' num2str(mdim) '_H.mat']);
save(goptofile);

% Cluster optimised hyperplanes

Lcluster(gopto,ctol,dopto);

% Save optimisation results


if histo
	optfile = fullfile(moddir,[optname '_mdim_' num2str(mdim) '_H.mat']);
else
	optfile = fullfile(moddir,[optname '_mdim_' num2str(mdim) '_N.mat']);
end
fprintf('*** saving optimisation results in ''%s''... ',optfile);
save(optfile);
fprintf('done\n\n');

% Plotting by calling python scripts

% gopto_str = mat2str(gopto); % Convert the MATLAB array to a string
% cmd = ['python', ' /Users/borjan/code/python/tvb/ssdigraphs.py', ' opt_plotting', ' --gopto ', gopto_str]
% subprocess.call(cmd)

% Call Python script to plot dopto
%system('python /Users/borjan/code/python/tvb/ssdigraphs.py plot_dopto');

% Getting localisation of macro:

nweight = zeros(n,size(Lopto,3));
for k = 1:size(Lopto,3)
	nweight(:,k) = 1-gmetricsx(Lopto(:,:,k));
end

% Plotting distance from causall graph co-ordinate system to hyperplane
% co-ordinate system

% k = 1; wgraph2dot(nweight(:,k),eweight,gfile,[],gvprog);

% Adding the plotting for gp_plucker: That looks like it might be working
% need to look into this propely

dist = 1 - gmetricsxx(Lopto(:,:,1));
%pluc_dist = 1 - plucker(Lopto(:,:,1), true);
% gptitle = sprintf('distances of locally-optimal %d-macro-subspace from %d-dim subspace in original system', mdim, mdim);
% gpstem = fullfile(tempdir, 'opt_dist');
% gp_plucker(dist, 9, mdim, gptitle, gpstem, gpterm, gpscale, gpfsize, gpplot)


if histo
	optfile = fullfile(moddir,[optname '_mdim_' num2str(mdim) '_H.mat']);
else
	optfile = fullfile(moddir,[optname '_mdim_' num2str(mdim) '_N.mat']);
end

fprintf('*** saving optimisation results in ''%s''... ',optfile);
save(optfile);
fprintf('done\n\n');

clear mdim
