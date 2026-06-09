% make axis vx and vkx
function vxyz = make_axis_spatial(Nx,sampling_distance,mode)
switch mode
    case '101'
        if mod(Nx,2) == 1
            vxyz = sampling_distance*linspace(-(Nx-1)/2,(Nx-1)/2,Nx); % Nx is odd, vector of spatial axis x
        else
            vxyz = sampling_distance*linspace(-Nx/2,Nx/2-1,Nx); % Nx is even
        end
    case '01'
        if mod(Nx,2) == 1
            vxyz = sampling_distance*linspace(0,Nx-1,Nx); % vector of spatial axis z
        else
            vxyz = sampling_distance*linspace(0,Nx-1,Nx); % vector of spatial axis z
        end
end