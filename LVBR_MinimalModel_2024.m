% LV corticothalamic model of perceptual awareness - visual rivalry
% Christopher Whyte

rng('shuffle')
close all; clear
set(0,'defaulttextinterpreter','latex')

%% Parameters

%----- simulation settings
T = 25000;         % length in ms
DT = .1;           % integration step in ms
sim_time = length(1:DT:T); % simulation length in steps of DT

%----- number of neurons
n_e = 90;     % excitatory LV cortical neurons
n_i = n_e;    % inhibitory LV cortical interneurons
n_th = 10;    % non-specific matrix thalamus neurons
n_total = n_e + n_i + n_th;

%----- connectivity

% E/I ring
e_e = 3.5; e_i = 1; i_e = 5;
sd_e_e = .5; sd_e_i = 2; sd_i_e = 2;
J_e_e = zeros(n_e,n_e); J_e_i = zeros(n_e,n_e); J_i_e = zeros(n_e,n_e);

% Cortical ring coupling
theta = linspace(0,2*pi,n_e);  % location on ring
for j = 1:n_e
    for k = 1:n_e
        % euclidean distance in coordinates of unit circle
        dist = sqrt(((cos(theta(j)) - cos(theta(k)))^2 + (sin(theta(j)) - sin(theta(k)))^2));
        J_e_e(j,k) = e_e*exp(-.5*(dist./sd_e_e).^2)./(sd_e_e*sqrt(2*pi));
        J_e_i(j,k) = e_i*exp(-.5*(dist./sd_e_i).^2)./(sd_e_i*sqrt(2*pi));
        J_i_e(j,k) = i_e*exp(-.5*(dist./sd_i_e).^2)./(sd_i_e*sqrt(2*pi));
    end
end

% ring -> th
e_th = 4;
J_e_th = zeros(n_th,n_e);
pocket_size = 9; counter = 1;
for j = 1:n_th
    J_e_th(j,counter:counter+pocket_size-1) = e_th;
    counter = counter + pocket_size;
end

% ring <- th
th_d = 10;
J_th_d = zeros(n_e,n_th);
pocket_size = 9; counter = 1;
for j = 1:n_th
    J_th_d(counter:counter+pocket_size-1,j) = th_d;
    counter = counter + pocket_size;
end

NMDA_E_E = 0.35; AMPA_E_E = 1.75;
% collect soma <-> soma connectivity into matrix 
J_exc_AMPA = zeros(n_total,n_total);
J_exc_NMDA = zeros(n_total,n_total);
J_inh = zeros(n_total,n_total);
% e <-> e (1:90,1:90)
J_exc_AMPA(1:n_e,1:n_e) = J_e_e*AMPA_E_E;
J_exc_NMDA(1:n_e,1:n_e) = J_e_e*NMDA_E_E;
% e -> i (91:180,1:90)
J_exc_AMPA(n_e+1:n_e+n_i,1:n_e) = J_e_i;
J_exc_NMDA(n_e+1:n_e+n_i,1:n_e) = J_e_i;
% e -> th (181:190,1:90)
J_exc_AMPA(n_e+n_i+1:n_total,1:n_e) = J_e_th;
% i -> e (1:90,91:180)
J_inh(1:n_e,n_e+1:n_e+n_i) = J_i_e;

% thalamus -> dendrite connectivity
J_dend = zeros(n_e,n_total);
J_dend(1:n_e,n_e+n_i+1:n_total) = J_th_d;

%----- neuron params

% -- General
% C = membane capicitance
% v_r = resting potential
% v_th = instantanious spike threshold
% v_peak = spike cut off
% k = scaling coefficient for membrane non-linearity
% a = adaptation time constant
% b = sensitivity to subthreshold oscillations
% c = voltage reset
% d = adaptation reset

% -- Dendrite specific
% pd = pulse amplitude of BAP
% ld = leak current coefficient

% L5b dendrite params
C_d = 170;
v_d_r = -70;
k_d = 1200;
p_d = 2600;
l_d = 170/7;
a_d = 1/30;
b_d = -13;
calcium_criterion = -30;

% dendrite <-> soma coupling 
BAP_delay = .5/DT;
BAP_width = 2/DT;

% L5b soma params
ve_peak = 50;
Ce = 150; ve_r = -75; ve_t = -45; ke = 2.5;
ae = 0.01; be = 5; ce = -65; de = 250;
burst_c_reset = -55; burst_d_reset = 150;

% basket cell params
vi_peak = 25;
Ci = 20; vi_r = -55; vi_t = -40; ki = 1;
ai = 0.15; bi = 8; ci = -55; di = 200;

% matrix thalamus params
vth_peak = 35;
Cth = 200; vth_r = -60; vth_t = -50; kth = 1.6;
ath = 0.01; bth = 15; cth = -60; dth = 10;
thal_inhbition_b = 0;
thal_voltage_switch = -65;
thal_idx = (n_e+n_i+1):n_total;

% collect params into vector
v_peak = [ve_peak.*ones(n_e,1); vi_peak.*ones(n_i,1); vth_peak.*ones(n_th,1)];
C = [Ce.*ones(n_e,1); Ci.*ones(n_i,1); Cth.*ones(n_th,1)];
v_r = [ve_r.*ones(n_e,1); vi_r.*ones(n_i,1); vth_r.*ones(n_th,1)];
v_t = [ve_t.*ones(n_e,1); vi_t.*ones(n_i,1); vth_t.*ones(n_th,1)];
k = [ke.*ones(n_e,1); ki.*ones(n_i,1); kth.*ones(n_th,1)];
a = [ae.*ones(n_e,1); ai.*ones(n_i,1); ath.*ones(n_th,1)];
b = [be.*ones(n_e,sim_time); bi.*ones(n_i,sim_time); bth.*ones(n_th,sim_time)];
c = [ce.*ones(n_e,sim_time); ci.*ones(n_i,sim_time); cth.*ones(n_th,sim_time)];
d = [de.*ones(n_e,sim_time); di.*ones(n_i,sim_time); dth.*ones(n_th,sim_time)];

%----- synapse params

% reverse potentials
E_exc = 0; E_inh = -75; E_adapt = -80;

% adaptation current
tau_decay_adapt = 2000; delta_adapt = 0.065;

% NMDA
tau_decay_NMDA = 100; 

% AMPA
tau_decay_AMPA = 6;

% GABAa
tau_decay_GABA = 6;

% coupling zone (decay const from: "Simulation of Postsynaptic
% Glutamate Receptors Reveals Critical Features of Glutamatergic
% Transmission")`    
tau_decay_coupling = 800;

%----- input into ring

% poisson process drive to soma
fr_L = 1350;                           % L firing rate in hertz
fr_R = 1350;                           % R firing rate in hertz

burnin_period = 500/DT;
burnin_L = linspace(0,fr_L,burnin_period);
burnin_R = linspace(0,fr_R,burnin_period);

% gaussian orientation selectivity
n = 1:n_e;
NL = 23; NR = NL+.5*n_e; sd = 18;
gaussian_selectivity_L = zeros(n_total,1);
gaussian_selectivity_R = zeros(n_total,1);
gaussian_selectivity_L(1:n_e,1) = exp(-((n-NL)/sd).^2)';
gaussian_selectivity_R(1:n_e,1) = exp(-((n-NR)/sd).^2)';

%% Model Lesions/Stimulation

Baclofen = zeros(n_e,1);
% Baclofen = Baclofen + 200;
% Baclofen(1:45) = Baclofen(1:45) + 200;
% Baclofen(1:45) = Baclofen(1:45) - 200;

thalamus_inhibition = zeros(n_total,1);
% thalamus_inhibition = thalamus_inhibition - 200;
% thalamus_inhibition(n_e+n_i+1:end,:) = thalamus_inhibition(n_e+n_i+1:end,:) - 200;

%% Initialise

%----- Initialise neuron
v_d = v_d_r.*ones(n_e,1); 
u_d = zeros(n_e,1);
v = v_r; 
u = zeros(n_total,1);

%----- Initialise synapse matrices
g_syn_NMDA = zeros(n_total,1);
g_syn_NMDA_d = zeros(n_e,1);
g_syn_AMPA = zeros(n_total,1);
g_syn_AMPA_d = zeros(n_e,1);
g_syn_GABA = zeros(n_total,1);
g_syn_adapt = zeros(n_total,1);
g_syn_coupling = zeros(n_e,1);

%----- Initialise external drive
I_ext_AMPA_L = 0; I_ext_AMPA_R = 0;
I_ext_NMDA_L = 0; I_ext_NMDA_R = 0;

%----- initialise spike storage containers
firings = []; BAP = zeros(n_e,sim_time);
v_d_store = nan(n_e,sim_time);
u_d_store = nan(n_e,sim_time);
v_store = nan(n_total,sim_time);
u_store = nan(n_total,sim_time);
I_adapt_store = nan(n_total,sim_time);

%% Simulation loop

disp('Simulation start'); tic
for t = 1:(sim_time-1)

    % display sim time
    if mod(t,1000/DT)==0
        disp(['t = ',num2str(t*DT)]);
    end

    %---------- find spikes
    fired = find(v>=v_peak);
    if ~isempty(fired)
        firings = [firings; t+0*fired, fired];
        v_store(fired,t-1) = v_peak(fired);
        % reset condition
        v(fired) = c(fired,t);
        % update adaptation variable
        u(fired) = u(fired) + d(fired,t);
        % ----- BAP
        L5b_soma_fired = fired(fired<=n_e);
        BAP(L5b_soma_fired,t+BAP_delay:t+BAP_delay+BAP_width) = 1;
        % ----- outgoing spikes
        % adaptation current
        g_syn_adapt(L5b_soma_fired) = g_syn_adapt(L5b_soma_fired) + delta_adapt;
        % -- e,i,th conductance
        % AMPA
        g_syn_AMPA = g_syn_AMPA + sum(J_exc_AMPA(:,fired),2);
        % NMDA
        g_syn_NMDA = g_syn_NMDA + sum(J_exc_NMDA(:,fired),2);
        % GABA
        g_syn_GABA = g_syn_GABA + sum(J_inh(:,fired),2);
        % -- th -> d conductance
        % AMPA
        g_syn_AMPA_d = g_syn_AMPA_d + sum(J_dend(:,fired),2);
        % NMDA
        g_syn_NMDA_d = g_syn_NMDA_d + sum(J_dend(:,fired),2);
        % thalamus dependent L5b apical<->soma coupling
        g_syn_coupling = g_syn_coupling + (1-g_syn_coupling).*sum(J_dend(:,fired)./th_d,2);
    end

    % dendrite storage 
    v_d_store(:,t) = v_d;
    u_d_store(:,t) = u_d;

    % soma  storage
    v_store(:,t) = v;
    u_store(:,t) = u;

    g_syn_coupling_store(:,t) = g_syn_coupling;
    g_syn_adapt_store(:,t) = g_syn_adapt;

    %---------- external drive 

    drive_L = fr_L;
    drive_R = fr_R;
    if t < burnin_period
       drive_L = burnin_L(t);
       drive_R = burnin_R(t);
    end 

    spikes_L = poissrnd(drive_L*(DT/1000),1);
    spikes_R = poissrnd(drive_R*(DT/1000),1);

    % add external spikes  
    I_ext_AMPA_L = I_ext_AMPA_L + spikes_L;
    I_ext_NMDA_L = I_ext_NMDA_L + spikes_L;
    I_ext_AMPA_R = I_ext_AMPA_R + spikes_R;
    I_ext_NMDA_R = I_ext_NMDA_R + spikes_R;

    % integrate external drive 
    I_ext_AMPA_L = I_ext_AMPA_L + DT.*(-I_ext_AMPA_L./tau_decay_AMPA);
    I_ext_NMDA_L = I_ext_NMDA_L + DT.*(-I_ext_NMDA_L./tau_decay_NMDA);
    I_ext_AMPA_R = I_ext_AMPA_R + DT.*(-I_ext_AMPA_R./tau_decay_AMPA);
    I_ext_NMDA_R = I_ext_NMDA_R + DT.*(-I_ext_NMDA_R./tau_decay_NMDA);

    % add external drive to each "eye" into single variable
    I_ext_AMPA = gaussian_selectivity_L.*I_ext_AMPA_L + gaussian_selectivity_R.*I_ext_AMPA_R;
    I_ext_NMDA = gaussian_selectivity_L.*I_ext_NMDA_L + gaussian_selectivity_R.*I_ext_NMDA_R;

    %---------- integrate synapses
    % AMPA
    g_syn_AMPA = g_syn_AMPA + DT.*(-g_syn_AMPA/tau_decay_AMPA);
    g_syn_AMPA_d = g_syn_AMPA_d + DT.*(-g_syn_AMPA_d/tau_decay_AMPA);
    % NMDA
    g_syn_NMDA = g_syn_NMDA + DT.*(-g_syn_NMDA/tau_decay_NMDA);
    g_syn_NMDA_d = g_syn_NMDA_d + DT.*(-g_syn_NMDA_d/tau_decay_NMDA);
    % GABA
    g_syn_GABA = g_syn_GABA + DT.*(-g_syn_GABA/tau_decay_GABA);
    % adaptation
    g_syn_adapt = g_syn_adapt + DT.*(-g_syn_adapt./tau_decay_adapt);
    % coupling zone
    g_syn_coupling = g_syn_coupling + DT.*(-g_syn_coupling./tau_decay_coupling);

    %---------- integrate L5b apical dendrites
    % conductance
    Mg_gate = (((v_d+80)./60).^2)./(1 + (((v_d+80)./60).^2));
    I_syn_dend = Mg_gate.*g_syn_NMDA_d.*(E_exc-v_d) + g_syn_AMPA_d.*(E_exc-v_d) + eps; % total condunctance
    % BAP probability dependent of couplin
    BAP_prob = rand(n_e,1) < g_syn_coupling;
    % membrane potential
    DV_d = (DT./C_d).*(-l_d*(v_d - v_d_r) + k_d.*F_dend(v_d) + BAP_prob.*p_d.*BAP(:,t) + u_d + I_syn_dend + Baclofen);
    DU_d = DT.*a_d.*(b_d.*(v_d - v_d_r) - u_d); 
    v_d = v_d + DV_d;
    u_d = u_d + DU_d;
    % if calcium occurs spike change reset conditions on soma of L5b
    calcium_threshold = v_d >= calcium_criterion;
    calcium_spike_prob = rand(n_e,1) < g_syn_coupling;
    calcium_spike = find((calcium_spike_prob == 1) & (calcium_threshold == 1));
    d(calcium_spike,t+1) = burst_d_reset;
    c(calcium_spike,t+1) = burst_c_reset;

    I_dend_store(:,t) = I_syn_dend + BAP_prob.*p_d.*BAP(:,t) + Baclofen;

    %---------- integrate L5b soma, Basket cell, Thalamus
    
    % conductance
    Mg_gate = (((v+80)./60).^2)./(1 + (((v+80)./60).^2));
    g_syn_tot = g_syn_AMPA + Mg_gate.*g_syn_NMDA + g_syn_GABA + g_syn_adapt + I_ext_AMPA + Mg_gate.*I_ext_NMDA + eps; % total condunctance
    E_tot = (g_syn_AMPA.*E_exc + Mg_gate.*g_syn_NMDA.*E_exc + g_syn_GABA.*E_inh + g_syn_adapt.*E_adapt + I_ext_AMPA.*E_exc + Mg_gate.*I_ext_NMDA.*E_exc)./g_syn_tot; % total reverse potential
    % adaptation properties of thalamus depend on voltage
    thal_voltage_idx = thal_idx(v(n_e+n_i+1:end) <= thal_voltage_switch);
    b(thal_voltage_idx,t) = thal_inhbition_b;
    % membrane potential 
    DV = (DT./C).*(k.*(v - v_r).*(v - v_t) - u + g_syn_tot.*E_tot + thalamus_inhibition);
    DU = DT.*a.*(b(:,t).*(v - v_r) - u);

    v = (v + DV)./(1 + g_syn_tot.*(DT./C));
    u = u + DU;

    g_syn_total(:,t) = g_syn_tot;
    E_total(:,t) = E_tot;

    I_adapt_store(:,t) = g_syn_adapt;

end
disp('Simulation end'); toc

%% Analysis

exc_neuron_idx = NL;
inh_neuron_idx = n_e + exc_neuron_idx;
thal_idx = n_e + n_i + round(exc_neuron_idx/9);

firings_e = firings(firings(:,2)<=n_e,:);
firings_i = firings(firings(:,2)>=n_e+1 & firings(:,2)<=n_e+n_i,:);
firings_th = firings(firings(:,2)>=n_e+n_i+1,:);

%---------- average firing rates

% define Gaussian smoother 
SmoothGaussWindow = @(x,n) conv(x,gausswin(n)./sum(gausswin(n)),'same');

% soma spike rate
firings = firings_e;

spikesLT = firings(1<= firings(:,2) & firings(:,2)<=45,1);
spikesRT = firings(45<firings(:,2) & firings(:,2)<=90,1);

binsize = 1/DT; % 1ms
tstep = 1:binsize:sim_time-1;

rateL_count = histcounts(spikesLT,tstep)*(1000/DT)/binsize/45; 
rateR_count = histcounts(spikesRT,tstep)*(1000/DT)/binsize/45;
rateL_smoothed = SmoothGaussWindow(rateL_count,250);
rateR_smoothed = SmoothGaussWindow(rateR_count,250);

error_bound = 5;
percieveL = rateL_smoothed > (rateR_smoothed + error_bound);
percieveR = rateR_smoothed > (rateL_smoothed + error_bound);

% smooth dist to criticality
LDistCrit = downsample(mean(I_dend_store(1:45,:)-538.911),1/DT);
RDistCrit = downsample(mean(I_dend_store(46:end,:)-538.911),1/DT);
LDistCrit_smoothed = SmoothGaussWindow(LDistCrit,250);
RDistCrit_smoothed = SmoothGaussWindow(RDistCrit,250);

% proportion past critical boundry
LDistCrit_boundry = downsample(mean((I_dend_store(1:45,:)-538.911) > 0),1/DT);
RDistCrit_boundry = downsample(mean((I_dend_store(46:end,:)-538.911) > 0),1/DT);

% smooth bursting regime
BuRegimeL = downsample(mean(d(1:45,:)==150),1/DT);
BuRegimeR = downsample(mean(d(46:90,:)==150),1/DT);
BuRegimeL_smoothed = SmoothGaussWindow(BuRegimeL,250);
BuRegimeR_smoothed = SmoothGaussWindow(BuRegimeR,250);

% smooth adaptation
I_adapt_L = downsample(mean(I_adapt_store(1:45,:)),1/DT);
I_adapt_R = downsample(mean(I_adapt_store(46:90,:)),1/DT);
I_adapt_L_smoothed = SmoothGaussWindow(I_adapt_L,250);
I_adapt_R_smoothed = SmoothGaussWindow(I_adapt_R,250);

% smooth adaptation regime
I_coupling_L = downsample(mean(g_syn_coupling_store(1:45,:)),1/DT);
I_coupling_R = downsample(mean(g_syn_coupling_store(46:90,:)),1/DT);
I_coupling_L_smoothed = SmoothGaussWindow(I_coupling_L,250);
I_coupling_R_smoothed = SmoothGaussWindow(I_coupling_R,250);

%---------- find switch points

% window
window = 300;
adapt_window = 1000;

% ---------- firing rate
rate_diff = rateL_smoothed - rateR_smoothed;
fAbove = rate_diff .* (rate_diff >= 0);
rate_switchpoints = find(diff(fAbove>0));
% remove last switch
rate_switchpoints = rate_switchpoints(3:length(rate_switchpoints)-2);
% grab time series around switch points
for idx = 1:length(rate_switchpoints)
    rateL_switchlocked(idx,:) = rateL_smoothed(rate_switchpoints(idx)-window:rate_switchpoints(idx)+window);
    rateR_switchlocked(idx,:) = rateR_smoothed(rate_switchpoints(idx)-window:rate_switchpoints(idx)+window);
end
% sort into L -> R and R - L
L2R = rateL_switchlocked(:,1)>=5;
% stack L2R and R2L
rate_sup_dom = mean([rateL_switchlocked(L2R,:); rateR_switchlocked(~L2R,:)]);
rate_dom_sup = mean([rateR_switchlocked(L2R,:); rateL_switchlocked(~L2R,:)]);

% ---------- dist to criticality
crit_diff = LDistCrit_smoothed - RDistCrit_smoothed;
fAbove = crit_diff .* (crit_diff >= 0);
crit_switchpoints = find(diff(fAbove>0));
% remove last switch
crit_switchpoints = crit_switchpoints(2:length(crit_switchpoints)-1);
% grab time series around rate switch points
for idx = 1:length(rate_switchpoints)
    critL_switchlocked(idx,:) = LDistCrit_smoothed(rate_switchpoints(idx)-window:rate_switchpoints(idx)+window);
    critR_switchlocked(idx,:) = RDistCrit_smoothed(rate_switchpoints(idx)-window:rate_switchpoints(idx)+window);
end
% stack L2R and R2L
crit_sup_dom = mean([critL_switchlocked(L2R,:); critR_switchlocked(~L2R,:)]);
crit_dom_sup = mean([critR_switchlocked(L2R,:); critL_switchlocked(~L2R,:)]);

% ---------- proportion bursting 
Bu_diff = BuRegimeL_smoothed - BuRegimeR_smoothed;
fAbove = Bu_diff .* (Bu_diff >= 0);
Bu_switchpoints = find(diff(fAbove>0));
Bu_switchpoints = Bu_switchpoints(2:length(Bu_switchpoints)-1);
% grab time series around rate switch points
for idx = 1:length(rate_switchpoints)
    BuL_switchlocked(idx,:) = BuRegimeL_smoothed(rate_switchpoints(idx)-window:rate_switchpoints(idx)+window);
    BuR_switchlocked(idx,:) = BuRegimeR_smoothed(rate_switchpoints(idx)-window:rate_switchpoints(idx)+window);
end
% stack L2R and R2L
Bu_sup_dom = mean([BuL_switchlocked(L2R,:); BuR_switchlocked(~L2R,:)]);
Bu_dom_sup = mean([BuR_switchlocked(L2R,:); BuL_switchlocked(~L2R,:)]);

% ---------- adaptation

for idx = 1:length(rate_switchpoints)
    I_adapt_L_switchlocked(idx,:) = I_adapt_L_smoothed(rate_switchpoints(idx)-adapt_window:rate_switchpoints(idx)+adapt_window);
    I_adapt_R_switchlocked(idx,:) = I_adapt_R_smoothed(rate_switchpoints(idx)-adapt_window:rate_switchpoints(idx)+adapt_window);
end

% stack L2R and R2L
Adapt_sup_dom = mean([I_adapt_L_switchlocked(L2R,:); I_adapt_R_switchlocked(~L2R,:)]);
Adapt_dom_sup = mean([I_adapt_R_switchlocked(L2R,:); I_adapt_L_switchlocked(~L2R,:)]);

% ---------- coupling
for idx = 1:length(rate_switchpoints)
    I_coupling_L_switchlocked(idx,:) = I_coupling_L_smoothed(rate_switchpoints(idx)-adapt_window:rate_switchpoints(idx)+adapt_window);
    I_coupling_R_switchlocked(idx,:) = I_coupling_R_smoothed(rate_switchpoints(idx)-adapt_window:rate_switchpoints(idx)+adapt_window);
end

% stack L2R and R2L
coupling_sup_dom = mean([I_coupling_L_switchlocked(L2R,:); I_coupling_R_switchlocked(~L2R,:)]);
coupling_dom_sup = mean([I_coupling_R_switchlocked(L2R,:); I_coupling_L_switchlocked(~L2R,:)]);

% ---------- compute stats of dominant population
Dom_Bu_L = BuRegimeL_smoothed(percieveL);
Dom_Bu_R = BuRegimeR_smoothed(percieveR);
Dom_Bu = [Dom_Bu_L, Dom_Bu_R];
mean_dom_Bu = mean(Dom_Bu);
std_dom_Bu = std(Dom_Bu);

Dom_crit_L = LDistCrit_smoothed(percieveL);
Dom_crit_R = RDistCrit_smoothed(percieveR);
Dom_crit = [Dom_crit_L, Dom_crit_R];
mean_dom_crit = mean(Dom_crit);
std_dom_crit = std(Dom_crit);

I_coupling_L = I_coupling_L_smoothed(percieveL);
I_coupling_R = I_coupling_R_smoothed(percieveR);
Dom_coupling = [I_coupling_L, I_coupling_R];
mean_Dom_coupling = mean(Dom_coupling);
std_Dom_coupling = std(Dom_coupling);

Dom_crit_L_boundry = LDistCrit_boundry(percieveL);
Dom_crit_R_boundry = RDistCrit_boundry(percieveR);
Dom_crit_boundry = [Dom_crit_L_boundry, Dom_crit_R_boundry];
mean_dom_crit_boundry = mean(Dom_crit_boundry);
std_dom_crit_boundry = std(Dom_crit_boundry);


%% Figures

time = 1:(sim_time); % time in seconds

figure(1)
hold on
plot(time/1000*DT,v_d_store(exc_neuron_idx,:),'b')
plot(time/1000*DT,v_store(exc_neuron_idx,:),'k')
ylabel('Neuron index')
ylabel('Time (s)')
ax = gca;
ax.FontSize = 25;
set(gca, 'FontName', 'Times')
title('L5 Soma')

figure(2)
plot(time/1000*DT,v_store(inh_neuron_idx,:),'k')
ax = gca;
ax.FontSize = 25;
set(gca, 'FontName', 'Times')
title('Basket Cell')

figure(3)
plot(time/1000*DT,v_store(thal_idx,:),'k')
ax = gca;
ax.FontSize = 25;
set(gca, 'FontName', 'Times')
title('Thalamus')

figure(4)
hold on
title('Raster: ttL5')
plot(firings_e(:,1)/10000,firings_e(:,2),'k.');
xlim([0,sim_time/10000])
ylim([0,n_e])
ylabel('Neuron index')
xlabel('Time (s)')
ax = gca;
ax.FontSize = 25;
ax.LineWidth = 1.5;
set(gca, 'FontName', 'Times')

figure(5)
hold on
title('Raster: L5 Basket cell')
plot(firings_i(:,1)/10000,firings_i(:,2),'k.');
xlim([0,sim_time/10000])
ylim([n_e,n_e+n_i])
ylabel('Neuron index')
xlabel('Time (s)')
ax = gca;
ax.FontSize = 25;
set(gca, 'FontName', 'Times')

figure(6)
hold on
title('Raster: Thalamus')
plot(firings_th(:,1)/10000,firings_th(:,2)-180,'k.');
xlim([0,sim_time/10000])
ylabel('Neuron index')
xlabel('Time (s)')
ax.LineWidth = 1.5;
ax = gca;
ax.FontSize = 15;
set(gca, 'FontName', 'Times')

time = -window:1:window;
adapt_time = -adapt_window:1:adapt_window;
yline = linspace(0,6,size(adapt_time,2));

f = figure(7);
subplot(3,1,1); hold on
title('Firing rate')
plot(time,rate_sup_dom,'Color', [101 200 150]/256,'linewidth',3)
plot(time,rate_dom_sup,'k','linewidth',3)
ylim([0,60])
ax = gca;
ax.FontSize = 15;
set(gca, 'FontName', 'Times')
box off
axis tight
subplot(3,1,2); hold on
title('Proportion bursting')
plot(time,Bu_sup_dom,'linewidth',3, 'Color', [80	180	202]/256)
plot(time,Bu_dom_sup,'k','linewidth',3)
ylim([0 .75])
ax = gca;
ax.FontSize = 15;
set(gca, 'FontName', 'Times')
box off
axis tight
ax.LineWidth = 1.5;
subplot(3,1,3); hold on
title('Distance to bifurcation')
plot(time,crit_sup_dom,'Color', [70 120 220]/256,'linewidth',3)
plot(time,crit_dom_sup,'k','linewidth',3)
ylim([-600,100])
ax = gca;
ax.FontSize = 15;
set(gca, 'FontName', 'Times')
box off
axis tight
ax.LineWidth = 1.5;

figure(8)
subplot(2,1,1); hold on
title('Coupling probability')
plot(adapt_time,coupling_sup_dom,'linewidth',3, 'Color', [60 60 256]/256)
plot(adapt_time,coupling_dom_sup,'k','linewidth',3)
ax = gca;
ax.FontSize = 15;
set(gca, 'FontName', 'Times')
box off
ylim([0 1.2])
subplot(2,1,2); hold on
title('Adaptation')
ax = gca;
ax.LineWidth = 1.5;
ax.FontSize = 15;
plot(adapt_time,yline, 'Color',[1 1 1])
plot(adapt_time,Adapt_sup_dom,'b','linewidth',3)
plot(adapt_time,Adapt_dom_sup,'k','linewidth',3)
ylim([0 6])
yticks([0 2 4 6])
set(gca, 'FontName', 'Times')
box off
axis tight


%% functions 

function x_out = F_dend(x_in)
    E_d = -38; D_d = 6;
    x_out = 1./(1 + exp((E_d - x_in)/D_d));
end 

