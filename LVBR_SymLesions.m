%% LV corticothalamic model of visual rivalry symmetric lesion simulations

% christopher Whyte 16/01/25

rng('shuffle')
close all
clear

%% simulation settings

simulation = 0;

% number of sims
num_sims =  30;
% length of simulation in ms
T = 30000; 
DT = .1;
sim_time = length(0:DT:T);  

% 0 = no stim, 1 = baclofen, 2 = opto, 3 = thal inhibition
lesion_type = [1 2 3];
% 0 = symmetric, 1 = asymmetric 
lesion_asymmetry = 0;

%% Get model params and stimulus

fr = 1350;
params = ThalamoCorticalParams(DT,sim_time);
lesion_strength_range = [0, 200, 400];

%% run simulations

if simulation == 1
    BR_SymLesion = {};
    counter = 0;
    for lesion_idx = 1:length(lesion_type)
        lesion = lesion_type(lesion_idx);
        for sim = 1:num_sims
            for lesion_strength_idx = 1:length(lesion_strength_range)
                counter = counter + 1;
                disp(['simulation: ',num2str(counter)]);
                lension_strength = lesion_strength_range(lesion_strength_idx);
                stimulus = BR_Stimulus(DT,params,fr,fr,lesion,0,lension_strength);
                [BR_SymLesion_firings{lesion_strength_idx,lesion_idx,sim} ...
                 BR_SymLesion_b{lesion_strength_idx,lesion_idx,sim} ...
                 BR_SymLesion_I_dend{lesion_strength_idx,lesion_idx,sim} ...
                 BR_SymLesion_I_adapt{lesion_strength_idx,lesion_idx,sim} ... 
                 BR_SymLesion_I_coupling{lesion_strength_idx,lesion_idx,sim}] = ThalamoCorticalNetSimulator(DT,sim_time,params,stimulus);
            end 
        end 
    end 
    save('BR_SymLesion_firings','BR_SymLesion_firings','-v7.3')
    save('BR_SymLesion_b','BR_SymLesion_b','-v7.3')
    save('BR_SymLesion_I_dend','BR_SymLesion_I_dend','-v7.3')
    save('BR_SymLesion_I_adapt','BR_SymLesion_I_adapt','-v7.3')
    save('BR_SymLesion_I_coupling','BR_SymLesion_I_coupling','-v7.3')
elseif simulation == 0
    load('BR_SymLesion_firings')
    load('BR_SymLesion_b')
    load('BR_SymLesion_I_dend')
    load('BR_SymLesion_I_adapt')
    load('BR_SymLesion_I_coupling')
end 


%% analysis of dominance durations

for lesion_idx = 1:length(lesion_type)
    for sim = 1:num_sims
        for lesion_strength_idx = 1:length(lesion_strength_range)
            firings = BR_SymLesion_firings{lesion_strength_idx,lesion_idx,sim};
            firings = firings(firings(:,2)<=params.n_e,:);
            [durationL{lesion_strength_idx,lesion_idx,sim},durationR{lesion_strength_idx,lesion_idx,sim},...
             durations{lesion_strength_idx,lesion_idx,sim},alternations(lesion_strength_idx,lesion_idx,sim)] = CalculateDominanceDurations(DT,firings,sim_time,params);
        end 
    end 
end 

alternations_mean = squeeze(mean(alternations*2,3));
alternations_err = squeeze(std(alternations*2,0,3)./sqrt(num_sims));

for lesion_idx = 1:length(lesion_type)
    for lesion_strength_idx = 1:length(lesion_strength_range)
        total_durations_LV4 = [];
        for sim = 1:num_sims
            total_durations_LV4 = [total_durations_LV4, durations{lesion_strength_idx,lesion_idx,sim}];
        end 
        total_durations_dist{lesion_strength_idx,lesion_idx} = total_durations_LV4;
        Levelt_durations_mean(lesion_strength_idx,lesion_idx) = mean(total_durations_LV4);
        Levelt_durations_err(lesion_strength_idx,lesion_idx) = std(total_durations_LV4)/sqrt(num_sims);
    end 
end 

%% analysis of dynamical features

burst_L = zeros(3,4,10,T-1); burst_R = zeros(3,4,10,T-1);
crit_L = zeros(3,4,10,T-1); crit_R = zeros(3,4,10,T-1);
adapt_L = zeros(3,4,10,T-1); adapt_R = zeros(3,4,10,T-1);
coup_L = zeros(3,4,10,T-1); coup_R = zeros(3,4,10,T-1);

for lesion_idx = 1:length(lesion_type)
    for lesion_strength_idx = 1:length(lesion_strength_range)
        for sim = 1:num_sims
            burst_L(lesion_idx,lesion_strength_idx,sim,:) = mean(BR_SymLesion_b{lesion_strength_idx,lesion_idx,sim}(1:45,1:T-1)==150,1);
            burst_R(lesion_idx,lesion_strength_idx,sim,:) = mean(BR_SymLesion_b{lesion_strength_idx,lesion_idx,sim}(46:end,1:T-1)==150,1);
            crit_L(lesion_idx,lesion_strength_idx,sim,:) = mean(BR_SymLesion_I_dend{lesion_strength_idx,lesion_idx,sim}(1:45,1:T-1) - 538.911,1);
            crit_R(lesion_idx,lesion_strength_idx,sim,:) = mean(BR_SymLesion_I_dend{lesion_strength_idx,lesion_idx,sim}(46:end,1:T-1) - 538.911,1);
            crit_L_prop(lesion_idx,lesion_strength_idx,sim,:) = mean((BR_SymLesion_I_dend{lesion_strength_idx,lesion_idx,sim}(1:45,1:T-1) - 538.911)>0,1);
            crit_R_prop(lesion_idx,lesion_strength_idx,sim,:) = mean((BR_SymLesion_I_dend{lesion_strength_idx,lesion_idx,sim}(46:end,1:T-1) - 538.911)>0,1);
            adapt_L(lesion_idx,lesion_strength_idx,sim,:) = mean(BR_SymLesion_I_adapt{lesion_strength_idx,lesion_idx,sim}(1:45,1:T-1),1);
            adapt_R(lesion_idx,lesion_strength_idx,sim,:) = mean(BR_SymLesion_I_adapt{lesion_strength_idx,lesion_idx,sim}(46:end,1:T-1),1);
            coup_L(lesion_idx,lesion_strength_idx,sim,:) = mean(BR_SymLesion_I_coupling{lesion_strength_idx,lesion_idx,sim}(1:45,1:T-1),1);
            coup_R(lesion_idx,lesion_strength_idx,sim,:) = mean(BR_SymLesion_I_coupling{lesion_strength_idx,lesion_idx,sim}(46:end,1:T-1),1);
        end 
    end 
end 

perceiveL = zeros(length(lesion_type),length(lesion_strength_range),num_sims,T-1);
perceiveR = zeros(length(lesion_type),length(lesion_strength_range),num_sims,T-1);

for lesion_idx = 1:length(lesion_type)
    for sim = 1:num_sims
        for lesion_strength_idx = 1:length(lesion_strength_range)
            firings = BR_SymLesion_firings{lesion_strength_idx,lesion_idx,sim};
            spikesLT = firings(1<= firings(:,2) & firings(:,2)<=45,1);
            spikesRT = firings(46<=firings(:,2) & firings(:,2)<=90,1);
            binsize = 1/.1; % 1ms
            tstep = 1:binsize:sim_time-1;
            rateL = histcounts(spikesLT,tstep)*(1000/DT)/binsize/45; 
            rateR = histcounts(spikesRT,tstep)*(1000/DT)/binsize/45;
            smoothBM = @(x,n) conv(x,gausswin(n)./sum(gausswin(n)),'same');
            rateL = smoothBM(rateL,250);
            rateR = smoothBM(rateR,250);
            error_bound = 5;
            perceiveL(lesion_idx,lesion_strength_idx,sim,:) = rateL > (rateR + error_bound);
            perceiveR(lesion_idx,lesion_strength_idx,sim,:) = rateR > (rateL + error_bound);
        end 
    end 
end 

dom_crit_L = zeros(length(lesion_type),length(lesion_strength_range),num_sims);
dom_crit_R = zeros(length(lesion_type),length(lesion_strength_range),num_sims);
dom_burst_L = zeros(length(lesion_type),length(lesion_strength_range),num_sims);
dom_burst_R = zeros(length(lesion_type),length(lesion_strength_range),num_sims);

% -- average distance to criticality and burst prob during dominance
for lesion_idx = 1:length(lesion_type)
    for lesion_strength_idx = 1:length(lesion_strength_range)
        for sim = 1:num_sims
            dom_idx_L = squeeze(perceiveL(lesion_idx,lesion_strength_idx,sim,:));
            dom_idx_R = squeeze(perceiveR(lesion_idx,lesion_strength_idx,sim,:));
            dom_crit_L(lesion_idx,lesion_strength_idx,sim) = mean(crit_L(lesion_idx,lesion_strength_idx,sim,dom_idx_L==1),4);
            dom_crit_R(lesion_idx,lesion_strength_idx,sim) = mean(crit_R(lesion_idx,lesion_strength_idx,sim,dom_idx_R==1),4);
            dom_crit_L_prop(lesion_idx,lesion_strength_idx,sim) = mean(crit_L_prop(lesion_idx,lesion_strength_idx,sim,dom_idx_L==1),4);
            dom_crit_R_prop(lesion_idx,lesion_strength_idx,sim) = mean(crit_R_prop(lesion_idx,lesion_strength_idx,sim,dom_idx_R==1),4);
            dom_burst_L(lesion_idx,lesion_strength_idx,sim) = mean(burst_L(lesion_idx,lesion_strength_idx,sim,dom_idx_L==1),4);
            dom_burst_R(lesion_idx,lesion_strength_idx,sim) = mean(burst_R(lesion_idx,lesion_strength_idx,sim,dom_idx_R==1),4);
            dom_coup_L(lesion_idx,lesion_strength_idx,sim) = mean(coup_L(lesion_idx,lesion_strength_idx,sim,dom_idx_L==1),4);
            dom_coup_R(lesion_idx,lesion_strength_idx,sim) = mean(coup_R(lesion_idx,lesion_strength_idx,sim,dom_idx_R==1),4);
        end 
    end 
end 

mean_crit = mean(cat(3,dom_crit_L,dom_crit_R),3);
std_crit = std(cat(3,dom_crit_L,dom_crit_R),0,3);

mean_crit_prop = mean(cat(3,dom_crit_L_prop,dom_crit_R_prop),3);
std_crit_prop = std(cat(3,dom_crit_L_prop,dom_crit_R_prop),0,3);

mean_burst = mean(cat(3,dom_burst_L,dom_burst_R),3);
std_burst = std(cat(3,dom_burst_L,dom_burst_R),0,3);

mean_coup = mean(cat(3,dom_coup_L,dom_coup_R),3);
std_coup = std(cat(3,dom_coup_L,dom_coup_R),0,3);

%% figures

f = figure(1); hold on
errorbar(lesion_strength_range,alternations_mean(:,2),alternations_err(:,2),'Color',[1 1 1], 'LineWidth',4)
errorbar(lesion_strength_range,alternations_mean(:,1),alternations_err(:,1),'Color',[202 61 131]/256, 'LineWidth',4)
yticks(0:10:50)
xticks([0 200 400])
ax = gca;
ax.FontSize = 25;
ax.LineWidth = 1.5;
set(gca, 'FontName', 'Times')
xlabel('Perturbation strength (pA)')
ylabel('Alternations (per min)')
axis padded
box off
f.Position = [822,345,392,420];

f = figure(2); hold on
ax = gca;
errorbar(lesion_strength_range,alternations_mean(:,2),alternations_err(:,2),'Color',[4 162 214]/256, 'LineWidth',4)
xticks([0 200 400])
ax.FontSize = 25;
ax.LineWidth = 1.5;
set(gca, 'FontName', 'Times')
xlabel('Perturbation strength (pA)')
ylabel('Alternations (per min)')
axis padded
box off
f.Position = [822,345,392,420];

f = figure(3); hold on
errorbar(lesion_strength_range,alternations_mean(:,2),alternations_err(:,2),'Color',[1 1 1], 'LineWidth',4)
errorbar(lesion_strength_range,alternations_mean(:,3),alternations_err(:,3),'Color',[181 96 53]/256, 'LineWidth',4)
yticks(0:10:50)
xticks([0 200 400])
ax = gca;
ax.FontSize = 25;
ax.LineWidth = 1.5;
set(gca, 'FontName', 'Times')
xlabel('Perturbation strength (pA)')
ylabel('Alternations (per min)')
axis padded
box off
f.Position = [822,345,392,420];


f = figure(4); hold on
ax = gca;
plot(lesion_strength_range,mean_crit(1,:),'Color',[202 61 131]/256, 'LineWidth',4)
plot(lesion_strength_range,mean_crit(2,:),'Color',[4 162 214]/256, 'LineWidth',4)
plot(lesion_strength_range,mean_crit(3,:),'Color',[181 96 53]/256, 'LineWidth',4)
xticks([0 200 400])
ax.FontSize = 25;
ax.LineWidth = 1.5;
set(gca, 'FontName', 'Times')
xlabel('Perturbation strength (pA)')
ylabel('Dist to $I_{B1}$ (pA)')
axis padded
box off
f.Position = [822,345,392,420];

f = figure(5); hold on
ax = gca;
plot(lesion_strength_range,mean_burst(1,:),'Color',[202 61 131]/256, 'LineWidth',4)
plot(lesion_strength_range,mean_burst(2,:),'Color',[4 162 214]/256, 'LineWidth',4)
plot(lesion_strength_range,mean_burst(3,:),'Color',[181 96 53]/256, 'LineWidth',4)
xticks([0 200 400])
ax.FontSize = 25;
ax.LineWidth = 1.5;
set(gca, 'FontName', 'Times')
xlabel('Perturbation strength (pA)')
ylabel('Proportion bursting')
axis padded
box off
f.Position = [822,345,392,420];

f = figure(6); hold on
ax = gca;
plot(lesion_strength_range,mean_crit_prop(1,:),'Color',[202 61 131]/256, 'LineWidth',4)
plot(lesion_strength_range,mean_crit_prop(2,:),'Color',[4 162 214]/256, 'LineWidth',4)
plot(lesion_strength_range,mean_crit_prop(3,:),'Color',[181 96 53]/256, 'LineWidth',4)
xticks([0 200 400])
ax.FontSize = 25;
ax.LineWidth = 1.5;
set(gca, 'FontName', 'Times')
xlabel('Perturbation strength (pA)')
ylabel('Proportion $> I_{B1}$ (pA)')
axis padded
box off
f.Position = [822,345,392,420];


f = figure(7);hold on
ax = gca;
plot(lesion_strength_range,Levelt_durations_mean(:,1)*sqrt(num_sims),'Color',[202 61 131]/256, 'linewidth',3)
plot(lesion_strength_range,Levelt_durations_mean(:,2)*sqrt(num_sims),'Color',[4 162 214]/256, 'linewidth', 3)
plot(lesion_strength_range,Levelt_durations_mean(:,3)*sqrt(num_sims),'Color',[181 96 53]/256, 'linewidth', 3)
ax.FontSize = 25;
ax.LineWidth = 1.5;
f.Position = [822,345,392,420];
set(gca, 'FontName', 'Times')
xlabel('Perturbation strength (pA)')
ylabel('Dominance duration $\sigma$')
axis padded
box off


f = figure(8);hold on
ax = gca;
plot(lesion_strength_range,Levelt_durations_err(:,1)*sqrt(num_sims),'Color',[202 61 131]/256, 'linewidth',3)
plot(lesion_strength_range,Levelt_durations_err(:,2)*sqrt(num_sims),'Color',[4 162 214]/256, 'linewidth', 3)
plot(lesion_strength_range,Levelt_durations_err(:,3)*sqrt(num_sims),'Color',[181 96 53]/256, 'linewidth', 3)
ax.FontSize = 25;
ax.LineWidth = 1.5;
f.Position = [822,345,392,420];
set(gca, 'FontName', 'Times')
xlabel('Perturbation strength (pA)')
ylabel('Dominance duration $\sigma$')
axis padded
box off

f = figure(9); hold on
ax = gca;
plot(lesion_strength_range,mean_coup(1,:),'Color',[202 61 131]/256, 'LineWidth',4)
plot(lesion_strength_range,mean_coup(2,:),'Color',[4 162 214]/256, 'LineWidth',4)
plot(lesion_strength_range,mean_coup(3,:),'Color',[181 96 53]/256, 'LineWidth',4)
xticks([0 200 400])
ax.FontSize = 25;
ax.LineWidth = 1.5;
set(gca, 'FontName', 'Times')
xlabel('Perturbation strength (pA)')
ylabel('Coupling propability (a.u.)')
axis padded
box off
f.Position = [822,345,392,420];

