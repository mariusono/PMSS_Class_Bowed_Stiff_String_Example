clc;
clear all;
% close all;

% Script controls
doNR = 1;
plotFlag = 0;
modulateBowPos = 1;
modulateOutPos = 1;

% General Params
% Fs = 44100/2;
Fs = 44100; 
k = 1/Fs;
dur = 4;
NF = floor(dur*Fs);
tVec = [0:1:NF-1].*k;

% String properties:
L = 0.5;
rho = 7850;
r = 5e-4;
T = 1000;
E = 2e11;
% sig0 = 0;
% sig1 = 0.005;

sig0 = 1;
% sig1 = 0.01;
sig1 = 0.005;

% sig0 = 0;
% sig1 = 0;

% Derived string params:
A = pi*(r)^2; % [m^2]
I = pi*r^4/4;
K = sqrt(E*I/(rho*A)); 

% % BOWING PARAMETERS STATIC MODEL:

% % Here we have big dif in doNR vs not doNR
% FB = 1.3; % dimensional case (it will be scaled with rho/A)
% vB = 0.1; % bow velocity (m/s)
% a = 80; % friction law free parameter (1/m^2) (a) % decreasing sig increases the stick time 

FB = 0.1; % [N]
% FB = 0.2; % [N]
vB = 0.2; % bow velocity (m/s)
a = 100; % friction law free parameter (1/m^2) (a) 
tol = 1e-7; % tolerance for Newton-Raphson method
A_NR = sqrt(2*a)*exp(1/2);

% input position
inp_bow_pos_x = 0.25; % in perc
bowPosVec = linspace(0.2,0.8,NF); % modulate position..

% output location
out_string_pos = 0.8;
outStringVec = linspace(0.8,0.2,NF);

% Target natural frequency
% N is number of intervals ! not points 
% f0 = 55; 
f0 = 220; 
% f0 = 440; 
% f0 = 880; 
c = f0*2*L; % wave speed term

% Stability condition
h_min = sqrt((c^2*k^2 + 4*sig1*k + sqrt((c^2*k^2 + 4*sig1*k)^2+16*K^2*k^2))/2);
N = floor(L/h_min); % no of intervals 
h = L/N;

% % Setup input force/velocity vector
durAttack = 0.2*dur;
durSustain = 0.5*dur;
durRelease = 0.05*dur;

linear_asr = generate_linear_asr_envelope_vector(durAttack,durSustain,durRelease,dur,Fs);
FB_vec = FB.*linear_asr;
vB_vec = vB.*linear_asr;

% Preallocate string displacement vectors
uPrev = zeros(N+1,1); % N+1 points and N intervals !
u = zeros(N+1,1); % same as uPrev but with some velocity ! 
uNext = zeros(N+1,1);

out = zeros(NF,1);
vel = zeros(NF,1);
disp = zeros(NF,1);
vrel_vec = zeros(NF,1);
count_vec = zeros(NF,1);
for n=1:NF

    FB = FB_vec(n);
    vB = vB_vec(n);

    if modulateBowPos
        bp = bowPosVec(n);
    else
        bp = inp_bow_pos_x;
    end

%     [I_B] = generate_interpolation_grid_1D(N+1,bp,'nearest'); % ratioed over the length..
%     [I_B] = generate_interpolation_grid_1D(N+1,bp,'linear'); % ratioed over the length..
    [I_B] = generate_interpolation_grid_1D(N+1,bp,'cubic'); % ratioed over the length..
    J_B = 1/h * I_B;

%     if mod(n,100) == 1
%         figure(123);
%         plot(I_B);
%     end

    if doNR
        I_B_J_B = sum(J_B(:).*J_B(:).* h);
    %     I_B_J_B = sum(I_B(:).*J_B(:)); % same ! 
    
        dt_min_u = (1/k).*(u-uPrev); 
        dxx_u = (1/h^2).*(u(3:N+1) - 2.*u(2:N) + u(1:N-1)); 
        dxxxx_u = (1/h^4).*(u(5:N+1) - 4.*u(4:N) + 6.*u(3:N-1) - 4.*u(2:N-2) + u(1:N-3)); 
        dt_min_dxx_u = (1/h^2).*(1/k).*(u(3:N+1) - 2.*u(2:N) + u(1:N-1) - (uPrev(3:N+1) - 2.*uPrev(2:N) + uPrev(1:N-1)));
    
        q = vB * (2/k + 2*sig0) -...
            (2/k) * sum(I_B(:).*dt_min_u(:)) - ... % i think u should use some circshifts here ? 
            (c^2) * sum(I_B(2:N).*dxx_u(:)) +...
            (K^2) * sum(I_B(3:N-1).*dxxxx_u(:)) -...
            (2*sig1) * sum(I_B(2:N).*dt_min_dxx_u(:)); 

        eps = 1;
        vrel_last = 0;
        count = 0;
        while eps>tol && count < 100
            vrel = vrel_last...
                   -((2/k + 2*sig0) * vrel_last ...
                   +(1/(rho*A))*(FB*I_B_J_B)*A_NR*vrel_last*exp(-a*vrel_last^2) + q)...
                   /((2/k + 2*sig0) + (1/(rho*A))*(FB*I_B_J_B)*A_NR*exp(-a*vrel_last^2)*(1+2*a*vrel_last^2));   
            eps = abs(vrel-vrel_last);
            vrel_last = vrel;
            count = count + 1;
            if count==99
                disp('NR fails');
            end
        end
        vrel = vrel_last;
        theta = A_NR.*vrel.*exp(-a.*vrel.^2);
        Fbow_vector = FB_vec(n) .* theta;       
        Fbow_vector = ones(size(u)).*Fbow_vector;
        
        vrel_vec(n) = vrel;
        count_vec(n) = count;
        Fbow_vec(n) = Fbow_vector(1);
    else
        vrel = (1/k) .* (u - uPrev) - vB_vec(n);
        vrel = sum(I_B(:).*vrel);
        theta = A_NR.*vrel.*exp(-a.*vrel.^2);
        Fbow_vector = FB_vec(n) .* theta;  
        vrel_vec(n) = vrel;   
        Fbow_vec(n) = Fbow_vector(1);        
    end

    uNext(3:N-1) = (k^2/(1+sig0*k))*( (c^2/h^2)*(u(4:N) - 2*u(3:N-1) + u(2:N-2))...
                    +(2*sig1/(k*h^2))*(u(4:N) - 2*u(3:N-1) + u(2:N-2) - uPrev(4:N) + 2*uPrev(3:N-1) - uPrev(2:N-2))...
                    - K^2*(1/h^4)*(u(5:N+1) - 4*u(4:N) + 6*u(3:N-1) - 4*u(2:N-2) + u(1:N-3))...
                    - J_B(3:N-1).*FB_vec(n) .* theta./(rho*A) ...
                    + (2/k^2)*u(3:N-1) - (1-sig0*k).*uPrev(3:N-1)./k^2  );

%     if strcmp(BC_string,'SimplySupported')
    uNext(2) = (k^2/(1+sig0*k))*( (c^2/h^2)*(u(3) - 2*u(2) + u(1))...
                    +(2*sig1/(k*h^2))*(u(3) - 2*u(2) + u(1) - uPrev(3) + 2*uPrev(2) - uPrev(1))...
                    - K^2*(1/h^4)*(u(4) - 4*u(3) + 6*u(2) - 4*u(1) - u(2))...
                    - J_B(2).*FB_vec(n) .* theta/(rho*A) ...
                    + (2/k^2)*u(2) - (1-sig0*k).*uPrev(2)./k^2  );                    

    uNext(N) = (k^2/(1+sig0*k))*( (c^2/h^2)*(u(N+1) - 2*u(N) + u(N-1))...
                    +(2*sig1/(k*h^2))*(u(N+1) - 2*u(N) + u(N-1) - uPrev(N+1) + 2*uPrev(N) - uPrev(N-1))...
                    - K^2*(1/h^4)*(-u(N) - 4*u(N+1) + 6*u(N) - 4*u(N-1) + u(N-2))...
                    - J_B(N).*FB_vec(n) .* theta./(rho*A) ...
                    + (2/k^2)*u(N) - (1-sig0*k).*uPrev(N)./k^2  );
%     end


    if modulateOutPos
        outPosPerc = outStringVec(n);
    else
        outPosPerc = out_string_pos;
    end

    % calculate velocity vector (state, or grid function) from u
    velVec = (1/(2*k) ) .* (uNext - uPrev);

    % cubic interpolation for the output to avoid clicks. check with
    % 'nearest' and see difference.
%     [I_out] = generate_interpolation_grid_1D(N+1,outPosPerc,'nearest'); % ratioed over the length..
    [I_out] = generate_interpolation_grid_1D(N+1,outPosPerc,'cubic'); % ratioed over the length..
    vel(n) = sum(I_out(:).*velVec(:));

%     vel(n) = (1/(2*k) ) .* (uNext(round(N * outPosPerc)) - uPrev(round(N * outPosPerc)));

    uPrev = u;
    u = uNext;

%     out(n) = u(round(N * outPosPerc));
    out(n) = sum(I_out(:).*u(:));

    
    if plotFlag
%         if n >= 16000
        if n >= 44000
%         if n >= 1

            vec = [1:length(u)];
            figure(123321123)
            hold off
            plot(vec.*h,u,'linewidth',2)
            grid on
            hold all
            plot(bp * L,sum(u.*I_B),'go','linewidth',3)
            plot(round(N * 0.8).*h,u(round(N * 0.8)),'ro','linewidth',2)  
%             ylim([-1,1].*10^-4)
%             pause
        end
    end
    
end

figure(1);
plot(tVec,out)
hold all


% Compare sound output as velocity vs displacement. Richer high freqs with
% vel output. Why ?
soundsc(vel,Fs);
% soundsc(out,Fs);


% % Experiment with outputing the sound at different sampling freq. That
% % stretches and pitches the sound. Play multiple sounds together.. Results
% % in some interesting stuff : ). Careful that I am not using soundsc here
% % so you should scale the sound to not blow up ur speakers. 
% sound(0.3.*vel./max(abs(vel)),Fs.*0.8);sound(0.4.*vel./max(abs(vel)),Fs.*0.5)
% sound(0.5.*vel./max(abs(vel)),Fs.*0.42);sound(0.5.*vel./max(abs(vel)),Fs.*0.56)
% 

% % Doing some freq analysis here.. You could try to do a spectrogram as well
% % and see what that shows. 
signal = vel./max(abs(vel));
% signal = out./max(abs(out));
% signal = vrel_vec./max(abs(vrel_vec));
% signal = signal./max(abs(signal));
% Nwin = length(signal)*2;
Nwin = length(signal);
xdft = fft(signal,Nwin);
xdft = xdft(1:Nwin/2+1);
psdx = (1/(Fs*Nwin)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/Nwin:Fs/2;

powerThing = 10*log10(psdx);

figure(767645);
hold all
plot(freq,powerThing,'linewidth',2)  % ref for units : https://stackoverflow.com/questions/1523814/units-of-a-fourier-transform-fft-when-doing-spectral-analysis-of-a-signal 
grid on
ylabel('Power/Frequency [dB/Hz]')
xlabel('Frequency [Hz]')
% set(gca,'xscale','log');


