% LV corticothalamic model of threshold detection and visual rivalry

% Christopher Whyte 10/01/25

function [firings, soma_burst_store, I_dend_store, I_adapt_store, I_coupling_store] = ThalamoCorticalNetSimulator(DT,sim_time,params,stimulus)
    
    % seed random number generator
    rng('shuffle');

    p = params;
    s = stimulus;
    
    clear params stimulus
    
    %---------- Initialise
    
    %----- Initialise neuron
    v_d = p.v_d_r.*ones(p.n_e,1); 
    u_d = zeros(p.n_e,1);
    v = p.v_r; 
    u = zeros(p.n_total,1);
    
    %----- Initialise synapse matrices
    g_syn_NMDA = zeros(p.n_total,1);
    g_syn_NMDA_d = zeros(p.n_e,1);
    g_syn_AMPA = zeros(p.n_total,1);
    g_syn_AMPA_d = zeros(p.n_e,1);
    g_syn_GABA = zeros(p.n_total,1);
    g_syn_adapt = zeros(p.n_total,1);
    g_syn_coupling = zeros(p.n_e,1);
    
    %----- Initialise external drive
    I_ext_AMPA_L = 0; I_ext_AMPA_R = 0;
    I_ext_NMDA_L = 0; I_ext_NMDA_R = 0;
    
    %----- initialise spike storage containers
    firings = []; BAP = zeros(p.n_e,sim_time-1);
    % v_d_store = zeros(p.n_e,sim_time-1);
    % u_d_store = zeros(p.n_e,sim_time-1);
    % v_store = zeros(p.n_total,sim_time-1);
    % u_store = zeros(p.n_total,sim_time-1);
    I_dend_store = zeros(p.n_e,sim_time-1);
    I_adapt_store = zeros(p.n_e,sim_time-1);
    I_coupling_store = zeros(p.n_e,sim_time-1);

    %---------- Simulation loop
    
    disp('Simulation start'); tic
    for t = 1:(sim_time-1)
    
        % display sim time
        if mod(t,1000/DT)==0
            disp(['t = ',num2str(t*DT)]);
        end
    
        %---------- find spikes
        fired = find(v>=p.v_peak);
        if ~isempty(fired)
            firings = [firings; t+0*fired, fired];
            % v_store(fired,t-1) = p.v_peak(fired);
            % reset condition
            v(fired) = p.c(fired,t);
            % update adaptation variable
            u(fired) = u(fired) + p.d(fired,t);
            % ----- BAP
            L5b_soma_fired = fired(fired<=p.n_e);
            BAP(L5b_soma_fired,t+p.BAP_delay:t+p.BAP_delay+p.BAP_width) = 1;
            % ----- outgoing spikes
            % adaptation current
            g_syn_adapt(L5b_soma_fired) = g_syn_adapt(L5b_soma_fired) + p.delta_adapt;
            % -- e,i,th conductance
            % AMPA
            g_syn_AMPA = g_syn_AMPA + sum(p.J_exc_AMPA(:,fired),2);
            % NMDA
            g_syn_NMDA = g_syn_NMDA + sum(p.J_exc_NMDA(:,fired),2);
            % GABA
            g_syn_GABA = g_syn_GABA + sum(p.J_inh(:,fired),2);
            % -- th -> d conductance
            % AMPA
            g_syn_AMPA_d = g_syn_AMPA_d + sum(p.J_dend(:,fired),2);
            % NMDA
            g_syn_NMDA_d = g_syn_NMDA_d + sum(p.J_dend(:,fired),2);
            % thalamus dependent L5b apical<->soma coupling
            g_syn_coupling = g_syn_coupling + (1-g_syn_coupling).*sum(p.J_dend(:,fired)./p.th_d,2);
        end
    
        % % dendrite storage 
        % v_d_store(:,t) = v_d;
        % u_d_store(:,t) = u_d;
        % 
        % % soma  storage
        % v_store(:,t) = v;
        % u_store(:,t) = u;
        % 
        % g_syn_coupling_store(:,t) = g_syn_coupling;
        % g_syn_adapt_store(:,t) = g_syn_adapt;

        %---------- external drive 

        drive_L = s.fr_L;
        drive_R = s.fr_R;
        if t < s.burnin_period
           drive_L = s.burnin_L(t);
           drive_R = s.burnin_R(t);
        end 

        spikes_L = poissrnd(drive_L*(DT/1000),1);
        spikes_R = poissrnd(drive_R*(DT/1000),1);

        % add external spikes  
        I_ext_AMPA_L = I_ext_AMPA_L + spikes_L;
        I_ext_NMDA_L = I_ext_NMDA_L + spikes_L;
        I_ext_AMPA_R = I_ext_AMPA_R + spikes_R;
        I_ext_NMDA_R = I_ext_NMDA_R + spikes_R;
    
        % integrate external drive 
        I_ext_AMPA_L = I_ext_AMPA_L + DT.*(-I_ext_AMPA_L./p.tau_decay_AMPA);
        I_ext_NMDA_L = I_ext_NMDA_L + DT.*(-I_ext_NMDA_L./p.tau_decay_NMDA);
        I_ext_AMPA_R = I_ext_AMPA_R + DT.*(-I_ext_AMPA_R./p.tau_decay_AMPA);
        I_ext_NMDA_R = I_ext_NMDA_R + DT.*(-I_ext_NMDA_R./p.tau_decay_NMDA);
    
        % add external drive to each "eye" into single variable
        I_ext_AMPA = s.gaussian_selectivity_L.*I_ext_AMPA_L + s.gaussian_selectivity_R.*I_ext_AMPA_R;
        I_ext_NMDA = s.gaussian_selectivity_L.*I_ext_NMDA_L + s.gaussian_selectivity_R.*I_ext_NMDA_R;
    
        %---------- integrate synapses
        % AMPA
        g_syn_AMPA = g_syn_AMPA + DT.*(-g_syn_AMPA/p.tau_decay_AMPA);
        g_syn_AMPA_d = g_syn_AMPA_d + DT.*(-g_syn_AMPA_d/p.tau_decay_AMPA);
        % NMDA
        g_syn_NMDA = g_syn_NMDA + DT.*(-g_syn_NMDA/p.tau_decay_NMDA);
        g_syn_NMDA_d = g_syn_NMDA_d + DT.*(-g_syn_NMDA_d/p.tau_decay_NMDA);
        % GABA
        g_syn_GABA = g_syn_GABA + DT.*(-g_syn_GABA/p.tau_decay_GABA);
        % adaptation
        g_syn_adapt = g_syn_adapt + DT.*(-g_syn_adapt./p.tau_decay_adapt);
        % coupling zone
        g_syn_coupling = g_syn_coupling + DT.*(-g_syn_coupling./p.tau_decay_coupling);
    
        %---------- integrate L5b apical dendrites
        % conductance
        Mg_gate = (((v_d+80)./60).^2)./(1 + (((v_d+80)./60).^2));
        I_syn_dend = Mg_gate.*g_syn_NMDA_d.*(p.E_exc-v_d) + g_syn_AMPA_d.*(p.E_exc-v_d) + eps; % total condunctance
        % BAP probability dependent on couplin
        BAP_prob = rand(p.n_e,1) < g_syn_coupling;
        % membrane potential
        DV_d = (DT./p.C_d).*(-p.l_d*(v_d - p.v_d_r) + p.k_d.*F_dend(v_d) + BAP_prob.*p.p_d.*BAP(:,t) + u_d + I_syn_dend + s.dend_perturbation);
        DU_d = DT.*p.a_d.*(p.b_d.*(v_d - p.v_d_r) - u_d); 
        v_d = v_d + DV_d;
        u_d = u_d + DU_d;
        % if calcium occurs spike change reset conditions on soma of L5b
        calcium_threshold = v_d >= p.calcium_criterion;
        calcium_spike_prob = rand(p.n_e,1) < g_syn_coupling;
        calcium_spike = find((calcium_spike_prob == 1) & (calcium_threshold == 1));
        p.d(calcium_spike,t+1) = p.burst_d_reset;
        p.c(calcium_spike,t+1) = p.burst_c_reset;

        I_dend_store(:,t) = I_syn_dend + s.dend_perturbation + BAP_prob.*p.p_d.*BAP(:,t);
        I_coupling_store(:,t) = g_syn_coupling;

        %---------- integrate L5b soma, Basket cell, Thalamus
        
        % conductance
        Mg_gate = (((v+80)./60).^2)./(1 + (((v+80)./60).^2));
        g_syn_tot = g_syn_AMPA + min(Mg_gate.*g_syn_NMDA,85) + g_syn_GABA + g_syn_adapt + I_ext_AMPA + min(Mg_gate.*I_ext_NMDA,85) + eps; % total condunctance
        E_tot = (g_syn_AMPA.*p.E_exc + min(Mg_gate.*g_syn_NMDA,85).*p.E_exc + g_syn_GABA.*p.E_inh + g_syn_adapt.*p.E_adapt + I_ext_AMPA.*p.E_exc + min(Mg_gate.*I_ext_NMDA,85).*p.E_exc)./g_syn_tot; % total reverse potential
        % adaptation properties of thalamus depend on voltage
        p.thal_voltage_idx = p.thal_idx(v(p.n_e+p.n_i+1:end) <= p.thal_voltage_switch);
        p.b(p.thal_voltage_idx,t) = p.thal_inhbition_b;
        % membrane potential 
        DV = (DT./p.C).*(p.k.*(v - p.v_r).*(v - p.v_t) - u + g_syn_tot.*E_tot + s.thal_perturbation);
        DU = DT.*p.a.*(p.b(:,t).*(v - p.v_r) - u);
    
        v = (v + DV)./(1 + g_syn_tot.*(DT./p.C));
        u = u + DU;

        I_adapt_store(:,t) = g_syn_adapt(1:p.n_e);
    
    end

    % down sample to 1ms
    I_dend_store = downsample(I_dend_store',1/DT)';
    I_adapt_store = downsample(I_adapt_store',1/DT)';
    I_coupling_store = downsample(I_coupling_store',1/DT)';
    soma_burst_store = downsample(p.d(1:p.n_e,:)',1/DT)';

    disp('Simulation end'); toc

end 

% Dendrite non-linearity 
function x_out = F_dend(x_in)
    E_d = -38; D_d = 6;
    x_out = 1./(1 + exp((E_d - x_in)/D_d));
end 
