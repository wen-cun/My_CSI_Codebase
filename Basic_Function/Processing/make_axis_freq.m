% make axis in frequency
function [vk] = make_axis_freq(N,sampling_distance,mode)
Fs = 1/sampling_distance;
F0 = Fs/N;
Kstep = F0;
Kmax = Fs/2;
switch mode
    case '101'
        if mod(N,2) == 1
            vk = linspace(-Kmax+Kstep/2,Kmax-Kstep/2,N); % vector of spatial freq axis -Kmax:0:+Kmax
        else
            vk = linspace(-Kmax,Kmax-Kstep,N);
        end
    case '01'
        if mod(N,2) == 1
            vk = linspace(0,Kmax-Kstep/2,round(N/2+.5)); % vector of spatial freq axis 0:+Kmax
        else
            vk = linspace(0,Kmax-Kstep,round(N/2));
        end
end