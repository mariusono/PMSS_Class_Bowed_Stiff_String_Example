function linear_asr = generate_linear_asr_envelope_vector(durAttack,durSustain,durRelease,durTotal,Fs)

ramp_fct = linspace(0,1,durAttack*Fs);
linear_asr_vec1 = interp1([0,1],[0,1],ramp_fct);

ramp_fct = linspace(0,1,durSustain*Fs);
linear_asr_vec2 = interp1([0,1],[1,1],ramp_fct);

ramp_fct = linspace(0,1,durRelease*Fs);
linear_asr_vec3 = interp1([0,1],[1,0],ramp_fct);

ramp_fct = linspace(0,1,durTotal*Fs-(durAttack*Fs+durSustain*Fs+durRelease*Fs));
linear_asr_vec4 = interp1([0,1],[0,0],ramp_fct);


linear_asr = [linear_asr_vec1,linear_asr_vec2,linear_asr_vec3,linear_asr_vec4];

if length(linear_asr)<durTotal*Fs
    linear_asr = [linear_asr,zeros(1,durTotal*Fs-length(linear_asr))];
end

end