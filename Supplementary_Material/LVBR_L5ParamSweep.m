%% LV corticothalamic rivalry "burstiness" parameter sweep simulation

% christopher Whyte 16/01/25

set(0,'defaulttextinterpreter','latex')
rng('shuffle')
close all
clear

%% simulation settings

simulation = 0;

% number of sims
num_sims = 10;
% length of simulation in ms
T = 30000; 
DT = .1;
sim_time = length(0:DT:T);  

%% Get model params and stimulus

fr = 1350; 
params = ThalamoCorticalParams(DT,sim_time);
stimulus = BR_Stimulus(DT,params,fr,fr,0,0,0); % 0 = no lesion

c_reset_range = [-52.5,-55, -57.5, -60, -62.5, -65];
d_reset_range = [125, 150, 175, 200, 225, 250];

%% run simulations

if simulation == 1
    counter = 0;
    for sim_idx = 1:num_sims
        for c_idx = 1:length(c_reset_range)
            for d_idx = 1:length(d_reset_range)
                params.burst_c_reset = c_reset_range(c_idx);
                params.burst_d_reset = d_reset_range(d_idx);
                counter = counter + 1;
                disp(['simulation: ',num2str(counter)]);
               [Levs_L5Sweep_firings{c_idx,d_idx,sim_idx},~,~,~,~] = ThalamoCorticalNetSimulator(DT,sim_time,params,stimulus);
            end
        end 
    end 
    save('Levs_L5Sweep_firings','Levs_L5Sweep_firings','-v7.3')
elseif simulation == 0
    load('Levs_L5Sweep_firings.mat')
end 

%% analysis

for sim_idx = 1:num_sims
    for c_idx = 1:length(c_reset_range)
        for d_idx = 1:length(d_reset_range)
            firings = Levs_L5Sweep_firings{c_idx,d_idx,sim_idx};
            firings = firings(firings(:,2)<=params.n_e,:);
            [durationL{c_idx,d_idx,sim_idx},durationR{c_idx,d_idx,sim_idx},...
            durations{c_idx,d_idx,sim_idx},alternations(c_idx,d_idx,sim_idx)] = CalculateDominanceDurations(DT,firings,sim_time,params);
        end 
    end 
end 

alternations_mean = mean(alternations,2);
alternations_err = std(alternations,0,2)./sqrt(num_sims);

for c_idx = 1:length(c_reset_range)
    for d_idx = 1:length(d_reset_range)
        total_durations = [];
        for sim_idx = 1:num_sims
            total_durations = [total_durations, durations{c_idx,d_idx,sim_idx}];
        end 
        total_durations_dist{c_idx,d_idx} = total_durations;
        Levelt_durations_mean(c_idx,d_idx) = mean(total_durations);
        Levelt_durations_err(c_idx,d_idx) = std(total_durations)/sqrt(num_sims);
    end 
end 

%% figures

[X,Y] = meshgrid(c_reset_range,d_reset_range);

f = figure(1); hold on
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 1.5;
contourf(X,Y,Levelt_durations_mean, "LevelList",0:500:3000)
set(gca, 'FontName', 'Times')
xlabel('Membrane potential reset (c)')
ylabel('After spike adaptation (d)')
set(gca, 'XDir','reverse')
a=colorbar;
a.Label.String = 'Dominance duration (ms)';
colormap('cool')
% f.Position = [822,345,280,420];


