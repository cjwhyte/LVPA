% LV corticothalamic model of threshold detection and visual rivalry model params

% Christopher Whyte 16/01/25

function [firings, soma_burst_store, I_dend_store, I_coupling_store] = ThalamoCorticalNetSimulatorTD(DT,sim_time,params,stimulus)
    
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
    I_ext_AMPA = zeros(p.n_total,1); I_ext_NMDA = zeros(p.n_total,1); 
    I_ext_AMPA_dend = zeros(p.n_e,1); I_ext_NMDA_dend = zeros(p.n_e,1);
    
    %----- initialise spike storage containers
    firings = []; BAP = zeros(p.n_e,sim_time-1);
    % v_d_store = zeros(p.n_e,sim_time);
    % u_d_store = zeros(p.n_e,sim_time);
    % v_store = zeros(p.n_total,sim_time);
    % u_store = zeros(p.n_total,sim_time);
    I_dend_store = zeros(p.n_e,sim_time-1);
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
        spikes_soma = poissrnd(s.background_drive_soma*(DT/1000),[p.n_e + p.n_i, 1]);
        spikes_dend = poissrnd(s.background_drive_dend*(DT/1000),[p.n_e, 1]);

        % add external spikes  
        I_ext_AMPA(1:p.n_e+p.n_i) = I_ext_AMPA(1:p.n_e+p.n_i) + spikes_soma;
        I_ext_NMDA(1:p.n_e+p.n_i) = I_ext_NMDA(1:p.n_e+p.n_i) + spikes_soma;
        I_ext_AMPA_dend = I_ext_AMPA_dend + spikes_dend;
        I_ext_NMDA_dend = I_ext_NMDA_dend + spikes_dend;
    
        % integrate external drive 
        I_ext_AMPA = I_ext_AMPA + DT.*(-I_ext_AMPA./p.tau_decay_AMPA);
        I_ext_NMDA = I_ext_NMDA + DT.*(-I_ext_NMDA./p.tau_decay_NMDA);
        I_ext_AMPA_dend = I_ext_AMPA_dend + DT.*(-I_ext_AMPA_dend./p.tau_decay_AMPA);
        I_ext_NMDA_dend = I_ext_NMDA_dend + DT.*(-I_ext_NMDA_dend./p.tau_decay_NMDA);
    
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
        I_syn_dend = Mg_gate.*g_syn_NMDA_d.*(p.E_exc-v_d) + g_syn_AMPA_d.*(p.E_exc-v_d) + Mg_gate.*I_ext_NMDA_dend.*(p.E_exc-v_d) + I_ext_AMPA_dend.*(p.E_exc-v_d) + eps; % total condunctance
        % BAP probability dependent of couplin
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

        I_dend_store(:,t) = I_syn_dend + BAP_prob.*p.p_d.*BAP(:,t) + s.dend_perturbation;
        I_coupling_store(:,t) = g_syn_coupling;
      
        %---------- integrate L5b soma, Basket cell, Thalamus
        
        % conductance
        Mg_gate = (((v+80)./60).^2)./(1 + (((v+80)./60).^2));
        g_syn_tot = g_syn_AMPA + Mg_gate.*g_syn_NMDA + g_syn_GABA + g_syn_adapt + s.W_FF.*I_ext_AMPA + s.W_FF.*Mg_gate.*I_ext_NMDA + eps; % total condunctance
        E_tot = (g_syn_AMPA.*p.E_exc + Mg_gate.*g_syn_NMDA.*p.E_exc + g_syn_GABA.*p.E_inh + g_syn_adapt.*p.E_adapt + s.W_FF.*I_ext_AMPA.*p.E_exc + s.W_FF.*Mg_gate.*I_ext_NMDA.*p.E_exc)./g_syn_tot; % total reverse potential
        % adaptation properties of thalamus depend on voltage
        p.thal_voltage_idx = p.thal_idx(v(p.n_e+p.n_i+1:end) <= p.thal_voltage_switch);
        p.b(p.thal_voltage_idx,t) = p.thal_inhbition_b;
        % membrane potential 
        DV = (DT./p.C).*(p.k.*(v - p.v_r).*(v - p.v_t) - u + g_syn_tot.*E_tot + s.amplitude.*s.I_ext_const(:,t) + s.thal_perturbation);
        DU = DT.*p.a.*(p.b(:,t).*(v - p.v_r) - u);
    
        v = (v + DV)./(1 + g_syn_tot.*(DT./p.C));
        u = u + DU;
    
        % g_syn_total_store(:,t) = g_syn_tot;
    
    end

    % down sample to 1ms
    I_dend_store = downsample(I_dend_store',1/DT)';
    soma_burst_store = downsample(p.d',1/DT)';
    I_coupling_store = downsample(I_coupling_store',1/DT)';

    disp('Simulation end'); toc


end 

% Dendrite non-linearity 
function x_out = F_dend(x_in)
    E_d = -38; D_d = 6;
    x_out = 1./(1 + exp((E_d - x_in)/D_d));
end 
