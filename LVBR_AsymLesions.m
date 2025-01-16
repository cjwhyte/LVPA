%% LV corticothalamic model of visual rivalry asymmetric lesion simulations

% christopher Whyte 16/01/25

rng('shuffle')
close all
clear
set(0,'defaulttextinterpreter','latex')

%% simulation settings

simulation = 0;

% number of sims
num_sims = 30;
% length of simulation in ms
T = 30000; 
DT = .1;
sim_time = length(0:DT:T);  

% 1 = baclofen, 2 = opto, 3 = thal inhibition
lesion_type = [1 2 3];

% 0 = symmetric, 1 = asymmetric 
lesion_asymmetry = 1;

%% Get model params and stimulus

fr = 1350;
params = ThalamoCorticalParams(DT,sim_time);
lesion_strength_range = [0, 200, 400];

%% run simulations

if simulation == 1
    counter = 0;
    for lesion_idx = 1:length(lesion_type)
        lesion = lesion_type(lesion_idx);
        for sim = 1:num_sims
            for lesion_strength_idx = 1:length(lesion_strength_range)
                counter = counter + 1;
                disp(['simulation: ',num2str(counter)]);
                lesion_strength = lesion_strength_range(lesion_strength_idx);
                stimulus = BR_Stimulus(DT,params,fr,fr,lesion,lesion_asymmetry,lesion_strength);
                [BR_AsymLesion_firings{lesion_strength_idx,lesion_idx,sim} ...
                 BR_AsymLesion_b{lesion_strength_idx,lesion_idx,sim} ...
                 BR_AsymLesion_I_dend{lesion_strength_idx,lesion_idx,sim} ...
                 BR_AsymLesion_I_adapt{lesion_strength_idx,lesion_idx,sim} ...
                 BR_AsymLesion_I_coupling{lesion_strength_idx,lesion_idx,sim}] = ThalamoCorticalNetSimulator(DT,sim_time,params,stimulus);
            end 
        end 
    end 
    save('BR_AsymLesion_firings','BR_AsymLesion_firings', '-v7.3')
    save('BR_AsymLesion_b','BR_AsymLesion_b', '-v7.3')
    save('BR_AsymLesion_I_dend','BR_AsymLesion_I_dend', '-v7.3')
    save('BR_AsymLesion_I_adapt','BR_AsymLesion_I_adapt', '-v7.3')
    save('BR_AsymLesion_I_coupling','BR_AsymLesion_I_coupling', '-v7.3')
elseif simulation == 0
    load('BR_AsymLesion_firings')
    load('BR_AsymLesion_b')
    load('BR_AsymLesion_I_dend')
    load('BR_AsymLesion_I_adapt')
    load('BR_AsymLesion_I_coupling')
end 


%% analysis of dominance durations

for lesion_idx = 1:length(lesion_type)
    for sim = 1:num_sims
        for lesion_strength_idx = 1:length(lesion_strength_range)
            firings = BR_AsymLesion_firings{lesion_strength_idx,lesion_idx,sim};
            firings = firings(firings(:,2)<=params.n_e,:);
            [durationL{lesion_strength_idx,lesion_idx,sim},durationR{lesion_strength_idx,lesion_idx,sim},...
             durations{lesion_strength_idx,lesion_idx,sim},alternations(lesion_strength_idx,lesion_idx,sim)] = CalculateDominanceDurations(DT,firings,sim_time,params);
        end 
    end 
end 

alternations_mean = mean(alternations,3);
alternations_err = std(alternations,0,3)./sqrt(num_sims);
for lesion_idx = 1:length(lesion_type)
    for lesion_strength_idx = 1:length(lesion_strength_range)
        durations_L = []; durations_R = []; 
        for sim = 1:num_sims
            durations_L = [durations_L, durationL{lesion_strength_idx,lesion_idx,sim}];
            durations_R = [durations_R, durationR{lesion_strength_idx,lesion_idx,sim}];
        end 
        total_duration_L{lesion_strength_idx,lesion_idx} = durations_L;
        total_duration_R{lesion_strength_idx,lesion_idx} = durations_R;
        total_durations_L_mean(lesion_strength_idx,lesion_idx) = mean(durations_L);
        total_durations_R_mean(lesion_strength_idx,lesion_idx) = mean(durations_R);
        total_durations_L_err(lesion_strength_idx,lesion_idx) = std(durations_L)./sqrt(num_sims);
        total_durations_R_err(lesion_strength_idx,lesion_idx) = std(durations_R)./sqrt(num_sims);
        total_durations_L_skew(lesion_strength_idx,lesion_idx) = kurtosis(durations_L);
        total_durations_R_skew(lesion_strength_idx,lesion_idx) = kurtosis(durations_R);

    end 
end 

no_lesion_duration_L = [total_duration_L{1,1},total_duration_L{1,2}, total_duration_L{1,3}];
no_lesion_duration_R = [total_duration_R{1,1},total_duration_R{1,2}, total_duration_R{1,3}];


%% analysis of dynamical features

burst_L = zeros(3,4,10,T-1); burst_R = zeros(3,4,10,T-1);
crit_L = zeros(3,4,10,T-1); crit_R = zeros(3,4,10,T-1);
crit_L_prop = zeros(3,4,10,T-1); crit_R_prop = zeros(3,4,10,T-1);
coup_L = zeros(3,4,10,T-1); coup_R = zeros(3,4,10,T-1);

for lesion_idx = 1:length(lesion_type)
    for lesion_strength_idx = 1:length(lesion_strength_range)
        for sim = 1:num_sims
            burst_L(lesion_idx,lesion_strength_idx,sim,:) = mean(BR_AsymLesion_b{lesion_strength_idx,lesion_idx,sim}(1:45,1:T-1)==150,1);
            burst_R(lesion_idx,lesion_strength_idx,sim,:) = mean(BR_AsymLesion_b{lesion_strength_idx,lesion_idx,sim}(46:90,1:T-1)==150,1);
            crit_L(lesion_idx,lesion_strength_idx,sim,:) = mean(BR_AsymLesion_I_dend{lesion_strength_idx,lesion_idx,sim}(1:45,1:T-1) - 538.911,1);
            crit_R(lesion_idx,lesion_strength_idx,sim,:) = mean(BR_AsymLesion_I_dend{lesion_strength_idx,lesion_idx,sim}(46:90,1:T-1) - 538.911,1);
            crit_L_prop(lesion_idx,lesion_strength_idx,sim,:) = mean((BR_AsymLesion_I_dend{lesion_strength_idx,lesion_idx,sim}(1:45,1:T-1) - 538.911)>0,1);
            crit_R_prop(lesion_idx,lesion_strength_idx,sim,:) = mean((BR_AsymLesion_I_dend{lesion_strength_idx,lesion_idx,sim}(46:90,1:T-1) - 538.911)>0,1);
            coup_L(lesion_idx,lesion_strength_idx,sim,:) = mean(BR_AsymLesion_I_coupling{lesion_strength_idx,lesion_idx,sim}(1:45,1:T-1),1);
            coup_R(lesion_idx,lesion_strength_idx,sim,:) = mean(BR_AsymLesion_I_coupling{lesion_strength_idx,lesion_idx,sim}(46:90,1:T-1),1);
        end 
    end 
end 

% clear variables to free up memory
clear BR_AsymLesion_b BR_AsymLesion_I_dend 

perceiveL = zeros(length(lesion_type),length(lesion_strength_range),num_sims,T-1);
perceiveR = zeros(length(lesion_type),length(lesion_strength_range),num_sims,T-1);

for lesion_idx = 1:length(lesion_type)
    for sim = 1:num_sims
        for lesion_strength_idx = 1:length(lesion_strength_range)
            firings = BR_AsymLesion_firings{lesion_strength_idx,lesion_idx,sim};
            spikesLT = firings(1<= firings(:,2) & firings(:,2)<=45,1);
            spikesRT = firings(46<=firings(:,2) & firings(:,2)<=90,1);
            binsize = 1/DT; % 1ms
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
dom_crit_L_prop = zeros(length(lesion_type),length(lesion_strength_range),num_sims);
dom_crit_R_prop = zeros(length(lesion_type),length(lesion_strength_range),num_sims);
dom_burst_L = zeros(length(lesion_type),length(lesion_strength_range),num_sims);
dom_burst_R = zeros(length(lesion_type),length(lesion_strength_range),num_sims);
dom_coup_L = zeros(length(lesion_type),length(lesion_strength_range),num_sims);
dom_coup_R = zeros(length(lesion_type),length(lesion_strength_range),num_sims);

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

mean_crit_L = mean(dom_crit_L,3);
mean_crit_R = mean(dom_crit_R,3);
std_crit_L = std(dom_crit_L,0,3);
std_crit_R = std(dom_crit_R,0,3);

mean_crit_L_prop = mean(dom_crit_L_prop,3);
mean_crit_R_prop = mean(dom_crit_R_prop,3);
std_crit_L_prop = std(dom_crit_L_prop,0,3);
std_crit_R_prop = std(dom_crit_R_prop,0,3);

mean_burst_L = mean(dom_burst_L,3);
mean_burst_R = mean(dom_burst_R,3);
std_burst_L = std(dom_burst_L,0,3);
std_burst_R = std(dom_burst_R,0,3);

mean_coup_L = mean(dom_coup_L,3);
mean_coup_R = mean(dom_coup_R,3);
std_coup_L = std(dom_coup_L,0,3);
std_coup_R = std(dom_coup_R,0,3);

%% Analysis of adaptation

adapt_L = zeros(length(lesion_strength_range),length(lesion_type),num_sims,T-1); 
adapt_R = zeros(length(lesion_strength_range),length(lesion_type),num_sims,T-1);

for lesion_idx = 1:length(lesion_type)
    for lesion_strength_idx = 1:length(lesion_strength_range)
        for sim = 1:num_sims
            adapt_L(lesion_strength_idx,lesion_idx,sim,:) = mean(BR_AsymLesion_I_adapt{lesion_strength_idx,lesion_idx,sim}(1:45,1:T-1),1);
            adapt_R(lesion_strength_idx,lesion_idx,sim,:) = mean(BR_AsymLesion_I_adapt{lesion_strength_idx,lesion_idx,sim}(46:90,1:T-1),1);
        end 
    end 
end

window = 100;

for sim = 1:num_sims
    for lesion_strength_idx = 1:length(lesion_strength_range)
        for lesion_idx = 1:length(lesion_type)
            firings = BR_AsymLesion_firings{lesion_strength_idx,lesion_idx,sim};
            spikesLT = firings(1<= firings(:,2) & firings(:,2)<=45,1);
            spikesRT = firings(46<=firings(:,2) & firings(:,2)<=90,1);
            binsize = 1/DT; % 1ms
            tstep = 1:binsize:sim_time-1;
            rateL = histcounts(spikesLT,tstep)*(1000/DT)/binsize/45; 
            rateR = histcounts(spikesRT,tstep)*(1000/DT)/binsize/45;
            smoothBM = @(x,n) conv(x,gausswin(n)./sum(gausswin(n)),'same');
            rateL = smoothBM(rateL,250);
            rateR = smoothBM(rateR,250);
            rate_diff = rateL - rateR;
            fAbove = rate_diff .* (rate_diff >= 0);
            rate_switchpoints = find(diff(fAbove>0));
            rate_switchpoints = rate_switchpoints(2:length(rate_switchpoints)-1);
            % grab time series around switch points
            for idx = 1:length(rate_switchpoints)
                rateL_switchlocked(idx,:) = rateL(rate_switchpoints(idx)-window:rate_switchpoints(idx)+window);
                adaptL_switchlocked(idx) = adapt_L(lesion_strength_idx,lesion_idx,sim,rate_switchpoints(idx));
                adaptR_switchlocked(idx) = adapt_R(lesion_strength_idx,lesion_idx,sim,rate_switchpoints(idx));
            end
            % sort into L -> R and R - L
            L2R = rateL_switchlocked(:,1)<5;
            % stack L2R and R2L
            adapt_sup_dom(lesion_strength_idx,lesion_idx,sim) = mean([adaptL_switchlocked(L2R), adaptR_switchlocked(~L2R)]);
            adapt_dom_sup(lesion_strength_idx,lesion_idx,sim) = mean([adaptR_switchlocked(L2R), adaptL_switchlocked(~L2R)]);
        end 
    end 
end 

mean_adapt_sup_dom = squeeze(mean(adapt_sup_dom,3));
mean_adapt_dom_sup = squeeze(mean(adapt_dom_sup,3));
std_adapt_sup_dom = squeeze(std(adapt_sup_dom,0,3))./sqrt(num_sims);
std_adapt_dom_sup = squeeze(std(adapt_dom_sup,0,3))./sqrt(num_sims);

%% figures

lesion_idx = 1;
f = figure(1); hold on
ax = gca;
% title('Baclofen')
errorbar(lesion_strength_range,total_durations_L_mean(:,lesion_idx),total_durations_L_err(:,lesion_idx), 'Color',[202 61 131]/256, 'LineWidth',4)
errorbar(lesion_strength_range,total_durations_R_mean(:,lesion_idx),total_durations_R_err(:,lesion_idx), 'Color',[202 61 131]/256, 'LineWidth',4, 'LineStyle','--')
ax.FontSize = 25;
ax.LineWidth = 1.5;
f.Position = [822,345,392,420];
set(gca, 'FontName', 'Times')
xlabel('Perturbation strength (pA)')
ylabel('Dominance duration')
axis padded
box off

lesion_idx = 2;
f = figure(2); hold on
ax = gca;
% title('Opto')
errorbar(lesion_strength_range, total_durations_L_mean(:,lesion_idx), total_durations_L_err(:,lesion_idx), 'Color',[4 162 214]/256, 'LineWidth',4)
errorbar(lesion_strength_range, total_durations_R_mean(:,lesion_idx), total_durations_R_err(:,lesion_idx), 'Color',[4 162 214]/256, 'LineWidth',4,'LineStyle','--')
ax.FontSize = 25;
ax.LineWidth = 1.5;
f.Position = [822,345,392,420];
set(gca, 'FontName', 'Times')
xlabel('Perturbation strength (pA)')
ylabel('Dominance duration')
axis padded
box off

lesion_idx = 3;
f = figure(3); hold on
% title('Thal inhibition')
ax = gca;
errorbar(lesion_strength_range, total_durations_L_mean(:,lesion_idx), total_durations_L_err(:,lesion_idx), 'Color',[181 96 53]/256, 'LineWidth',4)
errorbar(lesion_strength_range, total_durations_R_mean(:,lesion_idx), total_durations_R_err(:,lesion_idx), 'Color',[181 96 53]/256, 'LineWidth',4,'LineStyle','--')
xticks([0 200 400])
ax.FontSize = 25;
ax.LineWidth = 1.5;
f.Position = [822,345,392,420];
set(gca, 'FontName', 'Times')
xlabel('Perturbation strength (pA)')
ylabel('Dominance duration')
axis padded
box off

f = figure(4);hold on
ax = gca;
plot(lesion_strength_range,mean_crit_L(1,:),'Color',[202 61 131]/256, 'linewidth',3)
plot(lesion_strength_range,mean_crit_R(2,:),'Color',[4 162 214]/256, 'linewidth', 3,'Linestyle', '--')
plot(lesion_strength_range,mean_crit_L(3,:),'Color',[181 96 53]/256, 'linewidth', 3)
plot(lesion_strength_range,mean_crit_R(1,:),'Color',[202 61 131]/256, 'linewidth',3,'Linestyle', '--')
plot(lesion_strength_range,mean_crit_L(2,:),'Color',[4 162 214]/256, 'linewidth', 3)
plot(lesion_strength_range,mean_crit_R(3,:),'Color',[181 96 53]/256, 'linewidth', 3,'Linestyle', '--')
ax.FontSize = 25;
yticks([-400 -300 -200 -100 0 100 200 300])
ax.LineWidth = 1.5;
f.Position = [822,345,392,420];
set(gca, 'FontName', 'Times')
xlabel('Perturbation strength (pA)')
ylabel('Dist to $I_{B1}$ (pA)')
axis padded
box off

f = figure(5);hold on
ax = gca;
plot(lesion_strength_range,mean_crit_L_prop(1,:),'Color',[202 61 131]/256, 'linewidth',3)
plot(lesion_strength_range,mean_crit_R_prop(2,:),'Color',[4 162 214]/256, 'linewidth', 3,'Linestyle', '--')
plot(lesion_strength_range,mean_crit_L_prop(3,:),'Color',[181 96 53]/256, 'linewidth', 3)
plot(lesion_strength_range,mean_crit_R_prop(1,:),'Color',[202 61 131]/256, 'linewidth',3,'Linestyle', '--')
plot(lesion_strength_range,mean_crit_L_prop(2,:),'Color',[4 162 214]/256, 'linewidth', 3)
plot(lesion_strength_range,mean_crit_R_prop(3,:),'Color',[181 96 53]/256, 'linewidth', 3,'Linestyle', '--')
ax.FontSize = 25;
ax.LineWidth = 1.5;
f.Position = [822,345,392,420];
set(gca, 'FontName', 'Times')
xlabel('Perturbation strength (pA)')
ylabel('Proportion $>$ $I_{B1}$ (pA)')
axis padded
box off

f = figure(6);hold on
ax = gca;
plot(lesion_strength_range,mean_burst_L(1,:),'Color',[202 61 131]/256, 'linewidth',3)
plot(lesion_strength_range,mean_burst_R(2,:),'Color',[4 162 214]/256, 'linewidth', 3,'Linestyle', '--')
plot(lesion_strength_range,mean_burst_L(3,:),'Color',[181 96 53]/256, 'linewidth', 3)
plot(lesion_strength_range,mean_burst_R(1,:),'Color',[202 61 131]/256, 'linewidth',3,'Linestyle', '--')
plot(lesion_strength_range,mean_burst_L(2,:),'Color',[4 162 214]/256, 'linewidth', 3)
plot(lesion_strength_range,mean_burst_R(3,:),'Color',[181 96 53]/256, 'linewidth', 3,'Linestyle', '--')
ax.FontSize = 25;
ax.LineWidth = 1.5;
f.Position = [822,345,392,420];
set(gca, 'FontName', 'Times')
xlabel('Perturbation strength (pA)')
ylabel('Proportion bursting')
axis padded
box off

f = figure(7);hold on
ax = gca;
plot(lesion_strength_range,total_durations_L_err(:,1)*sqrt(num_sims),'Color',[202 61 131]/256, 'linewidth',3)
plot(lesion_strength_range,total_durations_L_err(:,2)*sqrt(num_sims),'Color',[4 162 214]/256, 'linewidth', 3)
plot(lesion_strength_range,total_durations_L_err(:,3)*sqrt(num_sims),'Color',[181 96 53]/256, 'linewidth', 3)
plot(lesion_strength_range,total_durations_R_err(:,1)*sqrt(num_sims),'Color',[202 61 131]/256, 'linewidth',3,'Linestyle', '--')
plot(lesion_strength_range,total_durations_R_err(:,2)*sqrt(num_sims),'Color',[4 162 214]/256, 'linewidth', 3,'Linestyle', '--')
plot(lesion_strength_range,total_durations_R_err(:,3)*sqrt(num_sims),'Color',[181 96 53]/256, 'linewidth', 3,'Linestyle', '--')
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
plot(lesion_strength_range,mean_coup_L(1,:),'Color',[202 61 131]/256, 'linewidth',3)
plot(lesion_strength_range,mean_coup_L(2,:),'Color',[4 162 214]/256, 'linewidth', 3,'Linestyle', '--')
plot(lesion_strength_range,mean_coup_L(3,:),'Color',[181 96 53]/256, 'linewidth', 3)
plot(lesion_strength_range,mean_coup_R(1,:),'Color',[202 61 131]/256, 'linewidth',3,'Linestyle', '--')
plot(lesion_strength_range,mean_coup_R(2,:),'Color',[4 162 214]/256, 'linewidth', 3)
plot(lesion_strength_range,mean_coup_R(3,:),'Color',[181 96 53]/256, 'linewidth', 3,'Linestyle', '--')
ax.FontSize = 25;
ax.LineWidth = 1.5;
f.Position = [822,345,392,420];
set(gca, 'FontName', 'Times')
xlabel('Perturbation strength (pA)')
ylabel('Coupling propability (a.u.)')
axis padded
box off
