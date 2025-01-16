%% Levelt's second law simulation scaled up model

% Christopher Whyte 16/01/25

set(0,'defaulttextinterpreter','latex')
rng('shuffle')
close all
clear

%% simulation settings

simulation = 0;

% number of sims
num_sims = 20;
% length of simulation in ms
T = 20000; 
DT = .1;
sim_time = length(0:DT:T);  

%% run simulations

stimulus_strength = flip(1400:100:1700);
fr_L = 1700;
if simulation == 1
    counter = 0;
    for sim = 1:num_sims
        params = ThalamoCorticalParams_ScaleDale(DT,sim_time);
        for stim_idx = 1:length(stimulus_strength)
            counter = counter + 1;
            disp(['simulation: ',num2str(counter)]);
            fr_R = stimulus_strength(stim_idx);
            stimulus = BR_Stimulus_ScaleDale(DT,params,fr_L,fr_R,0,0,0); % 0 = no lesion
            [Levs_SecondLaw_firings{stim_idx,sim}, ...
             ~,~,~,~] = ThalamoCorticalNetSimulator(DT,sim_time,params,stimulus);
        end 
    end 
    % save('Levs_SecondLaw_firings_ScaleDale','Levs_SecondLaw_firings','-v7.3')
elseif simulation == 0
    load('Levs_SecondLaw_firings_ScaleDale')
    params = ThalamoCorticalParams_ScaleDale(DT,sim_time);
    stimulus = BR_Stimulus(DT,params,fr_L,fr_L,0,0,0);
end 

%% analysis

for sim = 1:num_sims
    for stim_idx = 1:length(stimulus_strength)
        firings = Levs_SecondLaw_firings{stim_idx,sim};
        firings = firings(firings(:,2)<=params.n_e,:);
        [durationL{stim_idx,sim},durationR{stim_idx,sim},...
         durations{stim_idx,sim},alternations(stim_idx,sim)] = CalculateDominanceDurations(DT,firings,sim_time,params);
    end 
end 

alternations_mean = mean(alternations,1);
alternations_err = std(alternations,0,2)./sqrt(num_sims);

for stim_idx = 1:length(stimulus_strength)
    durations_L = []; durations_R = []; 
    for sim = 1:num_sims
        durations_L = [durations_L, durationL{stim_idx,sim}];
        durations_R = [durations_R, durationR{stim_idx,sim}];
    end 
    total_durations_L_mean(stim_idx) = mean(durations_L);
    total_durations_R_mean(stim_idx) = mean(durations_R);
    total_durations_L_err(stim_idx) = std(durations_L)./sqrt(num_sims);
    total_durations_R_err(stim_idx) = std(durations_R)./sqrt(num_sims);
end 

%% figures

stim_difference = stimulus.fr_L - stimulus_strength;

f = figure(1); hold on
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 1.5;
errorbar(stim_difference,total_durations_L_mean./1000,total_durations_L_err./1000, 'k-', 'LineWidth',4)
errorbar(stim_difference,total_durations_R_mean./1000,total_durations_R_err./1000, 'k--', 'LineWidth',4)
set(gca, 'FontName', 'Times')
xlabel('Stimulus difference (Hz)')
ylabel('Dominance Durations')
set(gca, 'XDir','reverse')
f.Position = [822,345,280,420];
axis padded

