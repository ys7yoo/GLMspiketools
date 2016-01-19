% Unit testing code to ensure that conditional intensity is computed the
% same way during simulation as during fitting

%% 1. Set parameters and display for a GLM % ==============================

dtStim = .01; % Bin size for stimulus (in seconds).  (Equiv to 100Hz frame rate)
dtSp = .001;  % Bin size for simulating model & computing likelihood (must evenly divide dtStim);
nkt = 30;    % Number of time bins in stimulus filter k
ttk = dtStim*(-nkt+1:0)';  % time relative to spike of stim filter taps
ggsim = makeSimStruct_GLM(nkt,dtStim,dtSp); % Create GLM structure with default params


%% 2. Simulate model and compute conditional intensity ====================
slen = 5000; % Stimulus length (frames);  More samples gives better fit
swid = 1;
Stim = randn(slen,swid);  %  Run model on long, binary stimulus

[tsp,vmem,ispk] = simGLM(ggsim,Stim);  % run model and return conditional intensity
nsp = length(tsp);

%% 3. Make param object with "true" params ================================
ggTrue = makeFittingStruct_GLM(dtStim,dtSp);
ggTrue.k = ggsim.k;   % insert stimulus filter
ggTrue.dc = ggsim.dc; % insert dc value
ggTrue.ihw = ggsim.ih; % insert spike-history filter
ggTrue.ihbas = speye(length(ggsim.ih));
ggTrue.tsp = tsp; % spike times
ggTrue.mask = [];

%% 4. Compute intensity again (during computation of log-likelihood) ======
% (if desired, compare rrT computed here with vmem computed above).

% To make it fail if you want to: 
% Stim(1) = Stim(1)+.1;

[logliTrue, rrT,tt,Iinj,Istm,Ih] = neglogli_GLM(ggTrue,Stim);

%% 5. Unit tests

% On total log-conditional intensity
assert(max(abs(vmem-Iinj))<1e-8,'Unit failed: log-conditional intensity consistency');

% On just the spike-history component of the conditional intensity
assert(max(abs(ispk-Ih))<1e-8,'Unit failed: spike-history filter output consistency');

fprintf('Unit passed: condIntensityConsistency\n');