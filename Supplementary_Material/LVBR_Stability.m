%% LV corticothalamic model of visual rivalry

% Christopher Whyte

rng('default')
close all; clear
set(0,'defaulttextinterpreter','latex')

%% Simulation type

% 1 = symmetric initial conditions
% 2 = asymmetric initial conditions
% 3 = symmetric initial conditions + perturbation
% 4 = aymmetric initial conditions + perturbation

simulation_type = 1;

%% Parameters

%----- simulation settings
T = 10000;         % length in ms
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

% collect soma <-> soma connectivity into matrix 
J_exc_AMPA = zeros(n_total,n_total);
J_exc_NMDA = zeros(n_total,n_total);
J_inh = zeros(n_total,n_total);
% e <-> e (1:90,1:90)
J_exc_AMPA(1:n_e,1:n_e) = J_e_e*1.75;
J_exc_NMDA(1:n_e,1:n_e) = J_e_e*0.35;
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
gamma = 1;  % coupling param 

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
NMDA_max = 85; % (nS)

% AMPA
tau_decay_AMPA = 6;

% GABAa
tau_decay_GABA = 6;

% coupling zone (decay const from: "Simulation of Postsynaptic
% Glutamate Receptors Reveals Critical Features of Glutamatergic
% Transmission")   
tau_decay_coupling = 800;

%----- input into ring
external_drive = 1400;

% gaussian orientation selectivity
n = 1:n_e;
NL = 23; NR = NL+.5*n_e; sd = 18;
gaussian_selectivity_L = zeros(n_total,1);
gaussian_selectivity_R = zeros(n_total,1);
gaussian_selectivity_L(1:n_e,1) = exp(-((n-NL)/sd).^2)';
gaussian_selectivity_R(1:n_e,1) = exp(-((n-NR)/sd).^2)';

I_ext = external_drive.*gaussian_selectivity_L + external_drive.*gaussian_selectivity_R;

pulse = zeros(n_total,sim_time);
if simulation_type >= 3
    pulse(NL,5000/DT:5000/DT+1/DT) = pulse(NL,5000/DT:5000/DT+1/DT) + 50;
end 


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
v = v_r.*ones(n_total,1); 
if simulation_type == 2 || simulation_type == 4
    v(1:45) = v_r(1:45) + 5;
end 
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
    % membrane potential
    DV_d = (DT./C_d).*(-l_d*(v_d - v_d_r) + k_d.*F_dend(v_d) + p_d.*BAP(:,t) + u_d + I_syn_dend + Baclofen);
    DU_d = DT.*a_d.*(b_d.*(v_d - v_d_r) - u_d); 
    v_d = v_d + DV_d;
    u_d = u_d + DU_d;
    % if calcium occurs spike change reset conditions on soma of L5b
    calcium_threshold = v_d >= calcium_criterion;
    calcium_spike = find(calcium_threshold == 1);
    d(calcium_spike,t+1) = burst_d_reset;
    c(calcium_spike,t+1) = burst_c_reset;

    I_dend_store(:,t) = I_syn_dend + p_d.*BAP(:,t) + Baclofen;

    %---------- integrate L5b soma, Basket cell, Thalamus
    
    % conductance
    Mg_gate = (((v+80)./60).^2)./(1 + (((v+80)./60).^2));
    g_syn_tot = g_syn_AMPA + min(Mg_gate.*g_syn_NMDA,NMDA_max) + g_syn_GABA + g_syn_adapt + eps; % total condunctance
    E_tot = (g_syn_AMPA.*E_exc + min(Mg_gate.*g_syn_NMDA,NMDA_max).*E_exc + g_syn_GABA.*E_inh + g_syn_adapt.*E_adapt)./g_syn_tot; % total reverse potential
    % adaptation properties of thalamus depend on voltage
    thal_voltage_idx = thal_idx(v(n_e+n_i+1:end) <= thal_voltage_switch);
    b(thal_voltage_idx,t) = thal_inhbition_b;
    % membrane potential 
    DV = (DT./C).*(k.*(v - v_r).*(v - v_t) - u + I_ext + pulse(:,t) + g_syn_tot.*E_tot + thalamus_inhibition);
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

%% Figures

time = (1:(size(rateL_smoothed,2)))*DT/100; % time in seconds

figure(1)
hold on
plot(time,rateL_smoothed,'Color', [101 200 150]/256,'linestyle','-','LineWidth',2)
plot(time,rateR_smoothed,'k', 'linestyle','-','LineWidth',2)
ylabel('Firing rate (Hz)')
xlabel('Time (s)')
ax = gca;
ax.LineWidth = 1.5;
ax.FontSize = 25;
set(gca, 'FontName', 'Times')

figure(2)
hold on
plot(rateL_smoothed,rateR_smoothed,'Color', [101 200 150]/256, 'linestyle','-','LineWidth',2)
ylabel('Firing rate R (Hz)')
xlabel('Firing rate L (Hz)')
xticks(0:10:50)
yticks(0:10:50)
ax = gca;
ax.LineWidth = 1.5;
ax.FontSize = 25;
set(gca, 'FontName', 'Times')

figure(3)
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

figure(4)
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

figure(5)
hold on
title('Raster: Thalamus')
plot(firings_th(:,1)/10000,firings_th(:,2)-180,'k.');
xlim([0,sim_time/10000])
ylabel('Neuron index')
xlabel('Time (s)')
ax.LineWidth = 1.5;
ax = gca;
ax.FontSize = 25;
set(gca, 'FontName', 'Times')



%% functions 

function x_out = F_dend(x_in)
    E_d = -38; D_d = 6;
    x_out = 1./(1 + exp((E_d - x_in)/D_d));
end 

