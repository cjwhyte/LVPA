% spike sorting script for LV thalamocortical model of binocularrivalry

function [durationL,durationR,Durations,Alternations] = CalculateDominanceDurations(DT,firings,sim_time,params)

    spikesLT = firings(1<= firings(:,2) & firings(:,2)<=(params.n_e/2),1);
    spikesRT = firings((params.n_e/2)+1<=firings(:,2) & firings(:,2)<=params.n_e,1);

    binsize = 1/DT; % 1ms
    tstep = 1:binsize:sim_time-1;
    
    rateL = histcounts(spikesLT,tstep)*(1000/DT)/binsize/(params.n_e/2); 
    rateR = histcounts(spikesRT,tstep)*(1000/DT)/binsize/(params.n_e/2);

    smoothBM = @(x,n) conv(x,gausswin(n)./sum(gausswin(n)),'same');
    rateL = smoothBM(rateL,250);
    rateR = smoothBM(rateR,250);

    error_bound = 5;
    perceiveL = rateL > (rateR + error_bound);
    perceiveR = rateR > (rateL + error_bound);

    vals = find(perceiveL>0);
    a = diff(vals);
    b = find([a inf]>1);
    durationL = diff([0 b]);

    vals = find(perceiveR>0);
    a = diff(vals);
    b = find([a inf]>1);
    durationR = diff([0 b]);

    Durations = [durationL(1:end-1), durationR(1:end-1)];
    Alternations = numel(durationL(1:end-1)) + numel(durationR(1:end-1));
   
end 
