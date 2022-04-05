function exp_curve_comb = generate_exponential_asr_envelope_vector(durAttack,durSustain,durRelease,durTotal,Fs)



% % Exponential attack-decay for soundsc
ramp_fct = linspace(0,1,durAttack*Fs);
exp_curve = exp(ramp_fct.*5)-1;
exp_curve = exp_curve./max(abs(exp_curve));

ramp_fct = linspace(0,1,durSustain*Fs);
exp_curve2 = interp1([0,1],[1,1],ramp_fct);

ramp_fct = linspace(0,1,durRelease*Fs);
exp_curve3 = exp(-ramp_fct.*3);
exp_curve3 = exp_curve3./max(abs(exp_curve3)).*(exp(-ramp_fct.*4)./max(abs(exp(-ramp_fct.*4))));

ramp_fct = linspace(0,1,durTotal*Fs-(durAttack*Fs+durSustain*Fs+durRelease*Fs));
exp_curve4 = interp1([0,1],[exp_curve3(end),0],ramp_fct);

exp_curve_comb = [exp_curve,exp_curve2,exp_curve3,exp_curve4];


if length(exp_curve_comb)<durTotal*Fs
    exp_curve_comb = [exp_curve_comb,zeros(1,durTotal*Fs-length(exp_curve_comb))];
end



end