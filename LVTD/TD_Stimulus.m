function stimulus = TD_Stimulus(DT,params,sim_time,lesion,lesion_magnitude)

    stimulus = {};

    % poisson process background drive
    stimulus.background_drive_soma = 600;                        
    stimulus.background_drive_dend = 50;    
    
    % amplitude of stimulus
    stimulus.amplitude = 1;

    % gaussian orientation selectivity
    N = 45; sd = 20; 
    stim_onset = 1000/DT; stim_width = 200/DT; 
    I_ext_const = zeros(params.n_total,sim_time); 
    for n = 1:params.n_e
        I_ext_const(n,stim_onset:stim_onset+stim_width) = exp(-((n-N)./sd).^2);
    end 
    stimulus.I_ext_const = I_ext_const;

    % FF weights
    stimulus.W_FF = ones(params.n_total,1);
    stimulus.W_FF(params.n_e+1:params.n_e+params.n_i,1) = 0.02.*stimulus.W_FF(params.n_e+1:params.n_e+params.n_i,1);

    dend_perturbation = zeros(params.n_e,1);
    thal_perturbation = zeros(params.n_total,1);
    if lesion == 1
        dend_perturbation = dend_perturbation - lesion_magnitude;
    elseif lesion == 2
        dend_perturbation = dend_perturbation + lesion_magnitude;
    elseif lesion == 3
       thal_perturbation(params.n_e+params.n_i+1:params.n_total) = thal_perturbation(params.n_e+params.n_i+1:params.n_total) - lesion_magnitude;
    end 

    stimulus.dend_perturbation = dend_perturbation;
    stimulus.thal_perturbation = thal_perturbation;

end 