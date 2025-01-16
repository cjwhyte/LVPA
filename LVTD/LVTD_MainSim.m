%% LV corticothalamic model threshold detection simulations

% christopher Whyte 16/01/24

rng('shuffle')
close all
clear

%% simulation settings

% set to 0 for preloading data and 1 to run sims
simulation = 1;

% number of sims
num_sims = 50;
% length of simulation in ms
T = 2000; 
DT = .1;
sim_time = length(0:DT:T);  

% 0 = no stim, 1 = baclofen, 2 = opto, 3 = thal inhibition
lesion_type = [0 1 2 3];

%% Get model params

params = ThalamoCorticalParams(DT,sim_time);

%% run simulations

pulse_amplitude = 0:50:600;
magnitude_range = 100:100:400;

if simulation == 1
    LVTD_firings = {}; LVTD_b = {};
    LVTD_I_dend = {}; LVTD_I_coupling = {};
    counter = 0;
    for magnitude_idx = 1:length(magnitude_range)
        lesion_magnitude = magnitude_range(magnitude_idx);
        for lesion_idx = 1:length(lesion_type)
            % get stimulus
            lesion = lesion_type(lesion_idx);
            stimulus = TD_Stimulus(DT,params,sim_time,lesion,lesion_magnitude);
            for sim = 1:num_sims
                for stim_idx = 1:length(pulse_amplitude)
                    counter = counter + 1;
                    disp(['simulation: ',num2str(counter)]);
                    stimulus.amplitude = pulse_amplitude(stim_idx);
                    [LVTD_firings{stim_idx,lesion_idx,magnitude_idx,sim} ...
                     LVTD_b{stim_idx,lesion_idx,magnitude_idx,sim} ...
                     LVTD_I_dend{stim_idx,lesion_idx,magnitude_idx,sim} ...
                     LVTD_I_coupling{stim_idx,lesion_idx,magnitude_idx,sim}] = ThalamoCorticalNetSimulatorTD(DT,sim_time,params,stimulus);
                end 
            end 
        end 
    end 
    save('LVTD_firings','LVTD_firings','-v7.3')
    save('LVTD_b','LVTD_b','-v7.3')
    save('LVTD_I_dend','LVTD_I_dend','-v7.3')
    save('LVTD_I_coupling','LVTD_I_coupling','-v7.3')
elseif simulation == 0
    load('LVTD_firings');
    load('LVTD_b');
    load('LVTD_I_dend');
    load('LVTD_I_coupling');
end 

%% Basic analysis and plots 

stim_idx = 3; lesion_idx = 1; magnitude_idx = 1; sim = 1; 

firings = LVTD_firings{stim_idx,lesion_idx,magnitude_idx,sim};

firings_e = firings(firings(:,2)<=params.n_e,:);
firings_i = firings(firings(:,2)>=params.n_e+1 & firings(:,2)<=params.n_e+params.n_i,:);
firings_th = firings(firings(:,2)>=params.n_e+params.n_i+1,:);

% ----- Average firing rate and CV for one sim

burn_in = 1000;

% soma
firings_exc = firings_e;
firings_inh = firings_i;

binsize = 1/DT; % 1ms
tstep = 1:binsize:sim_time-1;

rateE_count = histcounts(firings_exc,tstep)*(1000/DT)/binsize/90; 
rateI_count = histcounts(firings_inh,tstep)*(1000/DT)/binsize/90; 

rateTH_count = histcounts(firings_th,tstep)*(1000/DT)/binsize/10; 

SmoothGaussWindow = @(x,n) conv(x,gausswin(n)./sum(gausswin(n)),'same');

rateE_smoothed = SmoothGaussWindow(rateE_count,250);
rateI_smoothed = SmoothGaussWindow(rateI_count,250);

rateTH_smoothed = SmoothGaussWindow(rateTH_count,250);

time = linspace(0,2000,length(tstep)-1);

figure(1); hold on
plot(time(1,burn_in:end),rateE_smoothed(1,burn_in:end),'g','Linewidth',2)
plot(time(1,burn_in:end),rateI_smoothed(1,burn_in:end),'r','Linewidth',2)

isi = [];
for n = 1:params.n_e
    times_e = firings_e(firings_e(:,2) == n,1);
    times_i = firings_i(firings_i(:,2) == n,1);
    isi = [diff(times_e); isi,];
    isi = [diff(times_i); isi,];
    mean_isi_n(n) = mean(isi);
    std_isi_n(n) = std(isi);
end 

CV = mean(std_isi_n./mean_isi_n);

figure(2)
hold on
plot(firings_e(:,1)/(1000/DT),firings_e(:,2),'k.','MarkerSize',10);
xlim([0,sim_time/(1000/DT)])
ylim([0,params.n_e])
ylabel('Neuron index')
xlabel('Time (s)')
ax = gca;
ax.FontSize = 30; 
set(gca, 'FontName', 'Times')

%% Calculate psychometric and neurometric functions from spike counts

stim_onset = 1000/DT;
stim_offset = stim_onset + 1000/DT;
spike_counts = nan(num_sims,length(pulse_amplitude));
for magnitude_idx = 1:length(magnitude_range)
    for trl = 1:num_sims
        for intensity = 1:length(pulse_amplitude)
            % grab spikes: baseline
            firings = LVTD_firings{intensity,1,magnitude_idx,trl};
            firings = firings(firings(:,2)<=params.n_e,:);
            % grab spikes: baclofen
            firings_bac = LVTD_firings{intensity,2,magnitude_idx,trl};
            firings_bac = firings_bac(firings_bac(:,2)<=params.n_e,:);
            % grab spikes: opto
            firings_opto = LVTD_firings{intensity,3,magnitude_idx,trl};
            firings_opto = firings_opto(firings_opto(:,2)<=params.n_e,:);
            % grab spikes: thal
            firings_thal = LVTD_firings{intensity,4,magnitude_idx,trl};
            firings_thal = firings_thal(firings_thal(:,2)<=params.n_e,:);
            % post stimulus spike count: baseline
            spike_counts(trl,intensity,magnitude_idx) = size(firings(firings(:,1)>stim_onset & firings(:,1)<stim_offset,:),1);
            % post stimulus spike count: baclofen
            spike_counts_bac(trl,intensity,magnitude_idx) = size(firings_bac(firings_bac(:,1)>stim_onset & firings_bac(:,1)<stim_offset,1),1);
            % post stimulus spike count: opto
            spike_counts_opto(trl,intensity,magnitude_idx) = size(firings_opto(firings_opto(:,1)>stim_onset & firings_opto(:,1)<stim_offset,:),1);
            % post stimulus spike count: thal
            spike_counts_thal(trl,intensity,magnitude_idx) = size(firings_thal(firings_thal(:,1)>stim_onset & firings_thal(:,1)<stim_offset,:),1);
        end 
    end 
end 

% define criterion with respect to ~stimulation sims
min_spikes = min(spike_counts,[],"all");
max_spikes = max(spike_counts,[],"all");

criterion = min_spikes:1:max_spikes;
for crit = 1:length(criterion)
    noise_criterion(crit,:,:) = spike_counts(:,1,1) > criterion(crit);
    signal_criterion(crit,:,:,:) = spike_counts > criterion(crit);
    signal_criterion_bac(crit,:,:,:) = spike_counts_bac > criterion(crit);
    signal_criterion_opto(crit,:,:,:) = spike_counts_opto > criterion(crit);
    signal_criterion_thal(crit,:,:,:) = spike_counts_thal > criterion(crit);
end

noise = squeeze(mean(noise_criterion,2));
signal = squeeze(mean(signal_criterion,2));
signal_bac = squeeze(mean(signal_criterion_bac,2));
signal_opto = squeeze(mean(signal_criterion_opto,2));
signal_thal = squeeze(mean(signal_criterion_thal,2));

FA = noise;
H_avg = mean(signal(:,2:end),2);
miss = mean(signal(:,2:end)==0,2);

[~ ,threshold] = min(miss + FA);

% figure(); hold on
% plot(miss + FA)

% compute AUC
for magnitude_idx = 1:length(magnitude_range)
    for intensity = 1:length(pulse_amplitude)
        AUC(intensity,magnitude_idx) = trapz(flip(noise),flip(signal(:,intensity,magnitude_idx)));
        AUC_bac(intensity,magnitude_idx) = trapz(flip(noise),flip(signal_bac(:,intensity,magnitude_idx)));
        AUC_opto(intensity,magnitude_idx) = trapz(flip(noise),flip(signal_opto(:,intensity,magnitude_idx)));
        AUC_thal(intensity,magnitude_idx) = trapz(flip(noise),flip(signal_thal(:,intensity,magnitude_idx)));
    end 
end 

AUC = mean(AUC,2);

% compute psychometric curves
for magnitude_idx = 1:length(magnitude_range)
    for intensity = 1:length(pulse_amplitude)
        resp(:,intensity,magnitude_idx) = (spike_counts(:,intensity,magnitude_idx) > criterion(threshold));
        p_resp(:,intensity,magnitude_idx) = squeeze(mean(resp(:,intensity,magnitude_idx),1));
        % bac
        resp_bac(:,intensity,magnitude_idx) = (spike_counts_bac(:,intensity,magnitude_idx) > criterion(threshold));
        p_resp_bac(:,intensity,magnitude_idx) = squeeze(mean(resp_bac(:,intensity,magnitude_idx),1));
        % opto
        resp_opto(:,intensity,magnitude_idx) = (spike_counts_opto(:,intensity,magnitude_idx) > criterion(threshold));
        p_resp_opto(:,intensity,magnitude_idx) = squeeze(mean(resp_opto(:,intensity,magnitude_idx),1));
        % thal
        resp_thal(:,intensity,magnitude_idx) = (spike_counts_thal(:,intensity,magnitude_idx) > criterion(threshold));
        p_resp_thal(:,intensity,magnitude_idx) = squeeze(mean(resp_thal(:,intensity,magnitude_idx),1));
    end 
end 

p_resp = squeeze(mean(p_resp,3));

%% distance to criticality

crit_point = 538.911;
dist_crit = zeros(length(pulse_amplitude),length(lesion_type),length(magnitude_range),num_sims,round(T/2));
% average distance to criticality post stimulus period
for magnitude_idx = 1:length(magnitude_range)
    for lesion_idx = 1:length(lesion_type)
        for stim_idx = 1:length(pulse_amplitude)
            for sim = 1:num_sims
                dist_crit(stim_idx,lesion_idx,magnitude_idx,sim,:) = mean(LVTD_I_dend{stim_idx,lesion_idx,magnitude_idx,sim}(:,round(T/2)+1:end) - crit_point,1);
            end 
        end 
    end 
end 

tmean_dist_crit = squeeze(mean(dist_crit,5));
mean_dist_crit = squeeze(mean(tmean_dist_crit,4));
std_dist_crit = squeeze(std(tmean_dist_crit,0,4));

% average proportion bursting post stimulus period
Bu_Prop = zeros(length(pulse_amplitude),length(lesion_type),length(magnitude_range),num_sims,round(T/2));
% average distance to criticality post stimulus period
for magnitude_idx = 1:length(magnitude_range)
    for lesion_idx = 1:length(lesion_type)
        for stim_idx = 1:length(pulse_amplitude)
            for sim = 1:num_sims
                Bu_Prop(stim_idx,lesion_idx,magnitude_idx,sim,:) = mean(LVTD_b{stim_idx,lesion_idx,magnitude_idx,sim}(:,round(T/2)+1:end-1) == 150,1);
            end 
        end 
    end 
end 

tmean_Bu_Prop = squeeze(mean(Bu_Prop,5));
mean_Bu_Prop = squeeze(mean(tmean_Bu_Prop,4));
std_Bu_Prop = squeeze(std(tmean_Bu_Prop,0,4));

% average proportion bursting post stimulus period
I_coupling = zeros(length(pulse_amplitude),length(lesion_type),length(magnitude_range),num_sims,round(T/2));
% average distance to criticality post stimulus period
for magnitude_idx = 1:length(magnitude_range)
    for lesion_idx = 1:length(lesion_type)
        for stim_idx = 1:length(pulse_amplitude)
            for sim = 1:num_sims
                I_coupling(stim_idx,lesion_idx,magnitude_idx,sim,:) = mean(LVTD_I_coupling{stim_idx,lesion_idx,magnitude_idx,sim}(:,round(T/2)+1:end),1);
            end 
        end 
    end 
end 

tmean_I_coupling = squeeze(mean(I_coupling,5));
mean_I_coupling = squeeze(mean(tmean_I_coupling,4));
std_I_coupling = squeeze(std(tmean_I_coupling,0,4));

%% ----- plots

set(0,'defaulttextinterpreter','latex')
magnitude_idx = 3;

pulse_amplitude_plot = pulse_amplitude(:,1:end);

normalisation_const = min(AUC_thal(:,magnitude_idx));
normalisation_multiplier = max(AUC_opto(:,magnitude_idx)); 
proportionality_const = (normalisation_multiplier - normalisation_const);

yline = linspace(0,1,size(pulse_amplitude_plot,2));
xline = linspace(pulse_amplitude_plot(1),pulse_amplitude_plot (end),size(pulse_amplitude_plot ,2));

f = figure(7); hold on
axis tight
xticks([0 300 600])
yticks([0 .25 .5 .75 1])
ax = gca;
ax.FontSize = 25;
ax.LineWidth = 1.5;
set(gca,'TickDir','out'); 
plot(xline,yline, 'Color', [1 1 1])
plot(pulse_amplitude_plot,smooth((AUC_opto(:,magnitude_idx)-normalisation_const)*1/proportionality_const),'-', 'Color',[4 162 214]/256,'linewidth',3)
plot(pulse_amplitude_plot,smooth((AUC_opto(:,magnitude_idx)-normalisation_const)*1/proportionality_const),'.', 'Color',[4 162 214]/256,'linewidth',2, 'MarkerSize',40)
plot(pulse_amplitude_plot,smooth((AUC-normalisation_const)*1/proportionality_const),'-','Color',[117 113 113]/256,'linewidth',3)
plot(pulse_amplitude_plot,smooth((AUC-normalisation_const)*1/proportionality_const),'.','Color',[117 113 113]/256,'linewidth',2, 'MarkerSize',40)
set(gca, 'FontName', 'Times')
xlabel('Stim Intensity (pA)')
ylabel('P(response)')
f.Position = [822,345,392,420];
axis tight

f = figure(8); hold on
axis tight
xticks([0 300 600])
yticks([0 .25 .5 .75 1])
ax = gca;
ax.FontSize = 25;
ax.LineWidth = 1.5;
set(gca,'TickDir','out'); 
plot(xline,yline, 'Color', [1 1 1])
plot(pulse_amplitude_plot,smooth((AUC_bac(:,magnitude_idx)-normalisation_const)*1/proportionality_const),'-', 'Color', [202 61 131]/256,'linewidth',3)
plot(pulse_amplitude_plot,smooth((AUC_bac(:,magnitude_idx)-normalisation_const)*1/proportionality_const),'.', 'Color', [202 61 131]/256,'linewidth',2, 'MarkerSize',40)
plot(pulse_amplitude_plot,smooth((AUC-normalisation_const)*1/proportionality_const),'-','Color',[117 113 113]/256,'linewidth',3)
plot(pulse_amplitude_plot,smooth((AUC-normalisation_const)*1/proportionality_const),'.','Color',[117 113 113]/256,'linewidth',2, 'MarkerSize',40)
set(gca, 'FontName', 'Times')
xlabel('Stim Intensity (pA)')
ylabel('P(response)')
f.Position = [822,345,392,420];
axis tight

f = figure(9); hold on
axis tight
xticks([0 300 600])
yticks([0 .25 .5 .75 1])
ax = gca;
ax.FontSize = 25;
ax.LineWidth = 1.5;
set(gca,'TickDir','out'); 
plot(xline,yline, 'Color', [1 1 1])
plot(pulse_amplitude_plot,smooth((AUC_thal(:,magnitude_idx)-normalisation_const)*1/proportionality_const),'-', 'Color', [181 96 53]/256,'linewidth',3)
plot(pulse_amplitude_plot,smooth((AUC_thal(:,magnitude_idx)-normalisation_const)*1/proportionality_const),'.', 'Color', [181 96 53]/256,'linewidth',2, 'MarkerSize',40)
plot(pulse_amplitude_plot,smooth((AUC-normalisation_const)*1/proportionality_const),'-','Color',[117 113 113]/256,'linewidth',3)
plot(pulse_amplitude_plot,smooth((AUC-normalisation_const)*1/proportionality_const),'.','Color',[117 113 113]/256,'linewidth',2, 'MarkerSize',40)
set(gca, 'FontName', 'Times')
xlabel('Stim Intensity (pA)')
ylabel('P(response)')
f.Position = [822,345,392,420];
axis tight


f=figure(10); hold on
ax = gca;
ax.FontSize = 25;
ax.LineWidth = 1.5;
set(gca,'TickDir','out'); 
set(gca, 'FontName', 'Times')
xlabel('Stim Intensity (pA)')
ylabel('Spike count')
xticks([0 300 600])
plot(pulse_amplitude,mean(spike_counts(:,:,magnitude_idx),1)/params.n_e,'Color',[117 113 113]/256, 'linewidth',3)
plot(pulse_amplitude,mean(spike_counts_bac(:,:,magnitude_idx),1)/params.n_e,'Color',[202 61 131]/256, 'linewidth',3)
plot(pulse_amplitude,mean(spike_counts_opto(:,:,magnitude_idx),1)/params.n_e,'Color',[4 162 214]/256, 'linewidth',3)
plot(pulse_amplitude,mean(spike_counts_thal(:,:,magnitude_idx),1)/params.n_e,'Color',[181 96 53]/256, 'linewidth',3)
f.Position = [822,345,300,420];
axis tight

f=figure(11); hold on
ax = gca;
ax.FontSize = 25;
ax.LineWidth = 1.5;
set(gca,'TickDir','out'); 
set(gca, 'FontName', 'Times')
xlabel('Stim Intensity (pA)')
ylabel('Proportion bursting')
xticks([0 300 600])
plot(pulse_amplitude,squeeze(mean_Bu_Prop(:,1,magnitude_idx)),'Color',[117 113 113]/256, 'linewidth',3)
plot(pulse_amplitude,squeeze(mean_Bu_Prop(:,2,magnitude_idx)),'Color',[202 61 131]/256, 'linewidth',3)
plot(pulse_amplitude,squeeze(mean_Bu_Prop(:,3,magnitude_idx)),'Color',[4 162 214]/256, 'linewidth',3)
plot(pulse_amplitude,squeeze(mean_Bu_Prop(:,4,magnitude_idx)),'Color',[181 96 53]/256, 'linewidth',3)
f.Position = [822,345,300,420];
xlim([0,600])

f=figure(12); hold on
ax = gca;
ax.FontSize = 25;
ax.LineWidth = 1.5;
set(gca,'TickDir','out'); 
set(gca, 'FontName', 'Times')
xlabel('Stim Intensity (pA)')
ylabel('Dist to $I_{B1}$ (pA)')
xticks([0 300 600])
plot(pulse_amplitude,squeeze(mean_dist_crit(:,1,magnitude_idx)),'Color',[117 113 113]/256, 'linewidth',3)
plot(pulse_amplitude,squeeze(mean_dist_crit(:,2,magnitude_idx)),'Color',[202 61 131]/256, 'linewidth',3)
plot(pulse_amplitude,squeeze(mean_dist_crit(:,3,magnitude_idx)),'Color',[4 162 214]/256, 'linewidth',3)
plot(pulse_amplitude,squeeze(mean_dist_crit(:,4,magnitude_idx)),'Color',[181 96 53]/256, 'linewidth',3)
f.Position = [822,345,300,420];
xlim([0,600])

f=figure(13); hold on
ax = gca;
ax.FontSize = 25;
ax.LineWidth = 1.5;
set(gca,'TickDir','out'); 
set(gca, 'FontName', 'Times')
xlabel('Stim Intensity (pA)')
ylabel('Coupling propability (a.u.)')
xticks([0 300 600])
plot(pulse_amplitude,squeeze(mean_I_coupling(:,1,magnitude_idx)),'Color',[117 113 113]/256, 'linewidth',3)
plot(pulse_amplitude,squeeze(mean_I_coupling(:,2,magnitude_idx)),'Color',[202 61 131]/256, 'linewidth',3)
plot(pulse_amplitude,squeeze(mean_I_coupling(:,3,magnitude_idx)),'Color',[4 162 214]/256, 'linewidth',3)
plot(pulse_amplitude,squeeze(mean_I_coupling(:,4,magnitude_idx)),'Color',[181 96 53]/256, 'linewidth',3)
f.Position = [822,345,300,420];
ylim([.65,.85])
xlim([0,600])

% --- fit psychometric functions
x = (pulse_amplitude'+eps)/100;
stimulus_range = 0:.1:6;

logistic_intercept  = fittype("gamma + (1-gamma-lambda)/(1+exp(-beta*(x-alpha)))",...
                      dependent="y",independent="x",...
                      coefficients=["lambda" "beta" "alpha" "gamma"]);

f_control = fit(x,p_resp',logistic_intercept,'StartPoint',[0.003,2.5,2,0.01]); 
f_opto = fit(x,squeeze(p_resp_opto(:,:,magnitude_idx)'),logistic_intercept,'StartPoint',[0.003,2.5,2,0.01]);
f_bac = fit(x,squeeze(p_resp_bac(:,:,magnitude_idx)'),logistic_intercept,'StartPoint',[0.003,2.5,2,0.01]);
f_thal = fit(x,p_resp_thal(:,:,magnitude_idx)',logistic_intercept,'StartPoint',[0.003,2.5,2,0.01]);

% plot psychometric functions
f=figure(16); hold on
plot(x,p_resp,'.','MarkerSize',40,'Color',[117 113 113]/256)
plot(stimulus_range,f_control(stimulus_range),'Color',[117 113 113]/256, 'linewidth',3)
plot(x,p_resp_opto(:,:,magnitude_idx),'.','MarkerSize',40,'Color',[4 162 214]/256)
plot(stimulus_range,f_opto(stimulus_range),'Color',[4 162 214]/256, 'linewidth',3)
ylim([0,1.1])
ax = gca;
set(gca, 'FontName', 'Times')
set(gcf,'Color','w');
xlabel('Stimulus intensity (a.u.)')
ylabel('P(detection)')
ax.FontSize = 25;
ax.LineWidth = 1.5;
f.Position = [822,345,392,420];

f=figure(17); hold on
plot(x,p_resp,'.','MarkerSize',40,'Color',[117 113 113]/256)
plot(stimulus_range,f_control(stimulus_range),'Color',[117 113 113]/256, 'linewidth',3)
plot(x,p_resp_bac(:,:,magnitude_idx),'.','MarkerSize',40,'Color',[202 61 131]/256)
plot(stimulus_range,f_bac(stimulus_range),'Color',[202 61 131]/256, 'linewidth',3)
ylim([0,1.1])
ax = gca;
set(gca, 'FontName', 'Times')
set(gcf,'Color','w');
xlabel('Stimulus intensity (a.u.)')
ylabel('P(detection)')
ax.FontSize = 25;
ax.LineWidth = 1.5;
f.Position = [822,345,392,420];

f=figure(18); hold on
plot(x,p_resp,'.','MarkerSize',40,'Color',[117 113 113]/256)
plot(stimulus_range,f_control(stimulus_range),'Color',[117 113 113]/256, 'linewidth',3)
plot(x,p_resp_thal(:,:,magnitude_idx),'.','MarkerSize',40,'Color',[181 96 53]/256)
plot(stimulus_range,f_thal(stimulus_range),'Color',[181 96 53]/256, 'linewidth',3)
ylim([0,1.1])
ax = gca;
set(gca, 'FontName', 'Times')
set(gcf,'Color','w');
xlabel('Stimulus intensity (a.u.)')
ylabel('P(detection)')
ax.FontSize = 25;
ax.LineWidth = 1.5;
f.Position = [822,345,392,420];

