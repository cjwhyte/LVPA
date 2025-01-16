%% Levelt's fourth law simulation + distribution of dominance durations

% christopher Whyte 16/01/25

set(0,'defaulttextinterpreter','latex')
rng('shuffle')
close all
clear

%% simulation settings

simulation = 0;

% number of sims
num_sims = 40;
% length of simulation in ms
T = 30000; 
DT = .1;
sim_time = length(0:DT:T);  

%% Get model params and stimulus

params = ThalamoCorticalParams(DT,sim_time);

%% run simulations

stimulus_strength = 1300:50:1500;
if simulation == 1
    % -- Levelt's fourth law
    counter = 0;
    for sim = 1:num_sims
        for stim_idx = 1:length(stimulus_strength)
            counter = counter + 1;
            disp(['simulation: ',num2str(counter)]);
            fr_L = stimulus_strength(stim_idx);
            fr_R = stimulus_strength(stim_idx);
            stimulus = BR_Stimulus(DT,params,fr_L,fr_R,0,0,0); % 0 = no lesion
            [Levs_FourthLaw_firings{stim_idx,sim}, ~, ~, ~, ~] = ThalamoCorticalNetSimulator(DT,sim_time,params,stimulus);
        end 
    end 
    save('Levs_FourthLaw_firings','Levs_FourthLaw_firings','-v7.3')
elseif simulation == 0
    load('Levs_FourthLaw_firings')
end 


%% analysis

% -- Levelt's fourth law

for sim = 1:num_sims
    for stim_idx = 1:length(stimulus_strength)
        firings = Levs_FourthLaw_firings{stim_idx,sim};
        firings = firings(firings(:,2)<=params.n_e,:);
        [durationL{stim_idx,sim},durationR{stim_idx,sim},...
         durations{stim_idx,sim},alternations(stim_idx,sim)] = CalculateDominanceDurations(DT,firings,sim_time,params);
    end 
end 

alternations_mean = mean(alternations,2);
alternations_err = std(alternations,0,2)./sqrt(num_sims);

for stim_idx = 1:length(stimulus_strength)
    total_durations_LV4 = [];
    for sim = 1:num_sims
        total_durations_LV4 = [total_durations_LV4, durations{stim_idx,sim}];
    end 
    total_durations_dist{stim_idx} = total_durations_LV4;
    Levelt_durations_mean(stim_idx) = mean(total_durations_LV4);
    Levelt_durations_err(stim_idx) = std(total_durations_LV4)/sqrt(num_sims);
end 

%% distribution of dominance durations

total_durations_all = [];
for strength = 1:length(stimulus_strength)
    total_durations_all = [total_durations_all, total_durations_dist{strength}];
end 
total_durations = total_durations_all;
% exclude very short durations / fluctuations
total_durations = total_durations(total_durations>250);

dist_duration = fitdist(total_durations'./1000,'gamma'); % get params
dist_duration_norm = fitdist(total_durations'./1000,'normal'); % get params
dist_duration_lognorm = fitdist(total_durations'./1000,'lognormal'); % get params
LL1 = -negloglik(dist_duration);
LL2 = -negloglik(dist_duration_lognorm);
LL3 = -negloglik(dist_duration_norm);

%% figures

figure(1)
hold on
Nbins = histcounts(total_durations./1000');
h = histfit(total_durations./1000',length(Nbins), 'gamma');
y = get(gca, 'YTick'); % normalise y axis
set(gca, 'YTick', y, 'YTickLabel', round(y/numel(total_durations./1000'),2))
set(h(1),'facecolor',[76, 109, 172]/255); 
set(h(2),'color','k') % set colours
xlim([0,9])
ax = gca;
ax.FontSize = 25; 
set(gca, 'FontName', 'Times')
ylabel('Probability Density')
xlabel('Dominance Durations (s)')


f = figure(2);
ax = gca;
errorbar(stimulus_strength,alternations_mean,alternations_err,'k-', 'LineWidth',4)
ax.FontSize = 20;
ax.LineWidth = 1.5;
set(gca, 'FontName', 'Times')
xlabel('Stimulus strength (Hz)')
ylabel('Alternations (per min)')
f.Position = [822,345,280,420];
axis padded
box off
