function stimulus = BR_Stimulus(DT,params,fr_L,fr_R,lesion,lesion_asymmetry,lesion_strength)

    stimulus = {};

    % poisson process drive to soma (firing rate in hertz)
    stimulus.fr_L = fr_L;      
    stimulus.fr_R = fr_R;    

    stimulus.burnin_period = 500./DT;
    stimulus.burnin_L = linspace(0,stimulus.fr_L,stimulus.burnin_period);
    stimulus.burnin_R = linspace(0,stimulus.fr_R,stimulus.burnin_period);
   
    % gaussian orientation selectivity
    n = 1:params.n_e;
    stimulus.NL = 23; stimulus.NR = stimulus.NL + .5*params.n_e; stimulus.sd = 18;
    stimulus.gaussian_selectivity_L = zeros(params.n_total,1); 
    stimulus.gaussian_selectivity_R = zeros(params.n_total,1);
    stimulus.gaussian_selectivity_L(1:params.n_e,1) = exp(-((n-stimulus.NL)/stimulus.sd).^2)';
    stimulus.gaussian_selectivity_R(1:params.n_e,1) = exp(-((n-stimulus.NR)/stimulus.sd).^2)';

    % perturbation to model dendrites/thalamus
    dend_perturbation = zeros(params.n_e,1);
    thal_perturbation = zeros(params.n_total,1);
    if lesion_asymmetry == 0
        if lesion == 1
            dend_perturbation = dend_perturbation - lesion_strength;
        elseif lesion == 2
            dend_perturbation = dend_perturbation + lesion_strength;
        elseif lesion == 3
           thal_perturbation(params.n_e+params.n_i+1:params.n_total) = thal_perturbation(params.n_e+params.n_i+1:params.n_total) - lesion_strength;
        end 
    elseif lesion_asymmetry == 1
        if lesion == 1
            dend_perturbation(1:(params.n_e/2),1) = dend_perturbation(1:(params.n_e/2),1) - lesion_strength;
        elseif lesion == 2
            dend_perturbation(1:(params.n_e/2),1) = dend_perturbation(1:(params.n_e/2),1) + lesion_strength;
        elseif lesion == 3
           thal_perturbation((params.n_e + params.n_i + 1):(params.n_e + params.n_i + 1 + params.n_th/2),1)...
               = thal_perturbation((params.n_e + params.n_i + 1):(params.n_e + params.n_i + 1 + params.n_th/2),1) - lesion_strength;
        end 
    end 
    stimulus.dend_perturbation = dend_perturbation;
    stimulus.thal_perturbation = thal_perturbation;

end 