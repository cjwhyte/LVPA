function params = ThalamoCorticalParams(DT,sim_time)

    params = {};

    %----- number of neurons

    params.n_e = 90;     % excitatory LV cortical neurons
    params.n_i = params.n_e;    % inhibitory LV cortical interneurons
    params.n_th = 10;    % non-specific matrix thalamus neurons
    params.n_total = params.n_e + params.n_i + params.n_th;
    
    %----- connectivity
    
    % E/I ring
    params.e_e = 3.5; params.e_i = 1; params.i_e = 5;
    sd_e_e = .5; sd_e_i = 2; sd_i_e = 2;
    J_e_e = zeros(params.n_e,params.n_e); 
    J_e_i = zeros(params.n_e,params.n_e); 
    J_i_e = zeros(params.n_e,params.n_e);
    
    % Cortical ring coupling
    theta = linspace(0,2*pi,params.n_e);  % location on ring
    for j = 1:params.n_e
        for k = 1:params.n_e
            % euclidean distance in coordinates of unit circle
            dist = sqrt(((cos(theta(j)) - cos(theta(k)))^2 + (sin(theta(j)) - sin(theta(k)))^2));
            J_e_e(j,k) = params.e_e*exp(-.5*(dist./sd_e_e).^2)./(sd_e_e*sqrt(2*pi));
            J_e_i(j,k) = params.e_i*exp(-.5*(dist./sd_e_i).^2)./(sd_e_i*sqrt(2*pi));
            J_i_e(j,k) = params.i_e*exp(-.5*(dist./sd_i_e).^2)./(sd_i_e*sqrt(2*pi));
        end
    end
    
    % ring -> th
    params.e_th = 4;
    J_e_th = zeros(params.n_th,params.n_e);
    pocket_size = 9; counter = 1;
    for j = 1:params.n_th
        J_e_th(j,counter:counter+pocket_size-1) = params.e_th;
        counter = counter + pocket_size;
    end
    
    % ring <- th
    params.th_d = 10;
    J_th_d = zeros(params.n_e,params.n_th);
    pocket_size = 9; counter = 1;
    for j = 1:params.n_th
        J_th_d(counter:counter+pocket_size-1,j) = params.th_d;
        counter = counter + pocket_size;
    end
    
    % collect soma <-> soma connectivity into matrix 
    J_exc_AMPA = zeros(params.n_total,params.n_total);
    J_exc_NMDA = zeros(params.n_total,params.n_total);
    J_inh = zeros(params.n_total,params.n_total);
    % e <-> e (1:90,1:90)
    J_exc_AMPA(1:params.n_e,1:params.n_e) = J_e_e*1.75;
    J_exc_NMDA(1:params.n_e,1:params.n_e) = J_e_e*0.35;
    % e -> i (91:180,1:90)
    J_exc_AMPA(params.n_e+1:params.n_e+params.n_i,1:params.n_e) = J_e_i;
    J_exc_NMDA(params.n_e+1:params.n_e+params.n_i,1:params.n_e) = J_e_i;
    % e -> th (181:190,1:90)
    J_exc_AMPA(params.n_e+params.n_i+1:params.n_total,1:params.n_e) = J_e_th;
    % i -> e (1:90,91:180)
    J_inh(1:params.n_e,params.n_e+1:params.n_e+params.n_i) = J_i_e;
    
    % thalamus -> dendrite connectivity
    J_dend = zeros(params.n_e,params.n_total);
    J_dend(1:params.n_e,params.n_e+params.n_i+1:params.n_total) = J_th_d;
    
    params.J_exc_NMDA = J_exc_NMDA;
    params.J_exc_AMPA = J_exc_AMPA;
    params.J_inh = J_inh;
    params.J_dend = J_dend;
    
    %----- neuron params
    
    % -- General
    % C = membanecapicitance
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
    params.C_d = 170;
    params.v_d_r = -70;
    params.k_d = 1200;
    params.p_d = 2600;
    params.l_d = 170/7;
    params.a_d = 1/30;
    params.b_d = -13;
    params.calcium_criterion = -30;
    
    % dendrite <-> soma coupling 
    params.BAP_delay = .5/DT;
    params.BAP_width = 2/DT;
    
    % L5b soma params
    ve_peak = 50;
    Ce = 150; ve_r = -75; ve_t = -45; ke = 2.5;
    ae = 0.01; be = 5; ce = -65; de = 250;
    params.burst_c_reset = -55; params.burst_d_reset = 150;
    
    % basket cell params
    vi_peak = 25;
    Ci = 20; vi_r = -55; vi_t = -40; ki = 1;
    ai = 0.15; bi = 8; ci = -55; di = 200;
    
    % matrix thalamus params
    vth_peak = 35;
    Cth = 200; vth_r = -60; vth_t = -50; kth = 1.6;
    ath = 0.01; bth = 15; cth = -60; dth = 10;
    params.thal_inhbition_b = 0;
    params.thal_voltage_switch = -65;
    params.thal_idx = (params.n_e+params.n_i+1):params.n_total;
    
    % collect params into vector
    params.v_peak = [ve_peak.*ones(params.n_e,1); vi_peak.*ones(params.n_i,1); vth_peak.*ones(params.n_th,1)];
    params.C = [Ce.*ones(params.n_e,1); Ci.*ones(params.n_i,1); Cth.*ones(params.n_th,1)];
    params.v_r = [ve_r.*ones(params.n_e,1); vi_r.*ones(params.n_i,1); vth_r.*ones(params.n_th,1)];
    params.v_t = [ve_t.*ones(params.n_e,1); vi_t.*ones(params.n_i,1); vth_t.*ones(params.n_th,1)];
    params.k = [ke.*ones(params.n_e,1); ki.*ones(params.n_i,1); kth.*ones(params.n_th,1)];
    params.a = [ae.*ones(params.n_e,1); ai.*ones(params.n_i,1); ath.*ones(params.n_th,1)];
    params.b = [be.*ones(params.n_e,sim_time); bi.*ones(params.n_i,sim_time); bth.*ones(params.n_th,sim_time)];
    params.c = [ce.*ones(params.n_e,sim_time); ci.*ones(params.n_i,sim_time); cth.*ones(params.n_th,sim_time)];
    params.d = [de.*ones(params.n_e,sim_time); di.*ones(params.n_i,sim_time); dth.*ones(params.n_th,sim_time)];
    
    %----- synapse params
    
    % reverse potentials
    params.E_exc = 0; params.E_inh = -75; params.E_adapt = -80;
 
    % adaptation current
    params.tau_decay_adapt = 2000; params.delta_adapt = 0.065;
    
    % NMDA
    params.tau_decay_NMDA = 100;
    
    % AMPA
    params.tau_decay_AMPA = 6;
    
    % GABAa
    params.tau_decay_GABA = 6;
    
    % coupling zone (decay const from: "Simulation of Postsynaptic
    % Glutamate Receptors Reveals Critical Features of Glutamatergic
    % Transmission")`    
    params.tau_decay_coupling = 800;

end 