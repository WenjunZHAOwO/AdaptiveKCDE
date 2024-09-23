function [hx, N_effective_z,varxz,sigmaxz] = rule_of_thumb_x(zstar, hz, x,z,Gz,hzstar,k)

D = size(x,1);
N = size(x,2);


if nargin < 6
    hzstar = hz;
end

if nargin < 7
    k = 2;
end


 

% Kernel matrix in z
Kz =  Gz(zstar,hz,z);


Kzfull = Gz(z, hz, z);


% Number of effective neighbors in z space
% size(zstar)
% size(Kz)
N_effective_z =  sum( Kz./diag(Gz(zstar(1),hzstar, zstar(1))) );

% Estimate conditional mean and SD given z using kernel regression

xbarz =   x*( Kzfull./sum(Kzfull,1));

varxz =  ( x-xbarz).^2*( Kz./sum(Kz,1));


sigmaxz =  sqrt(varxz );

hx =  (4/(D+2))^(1/(D+4))*sigmaxz.* N_effective_z.^(-1/(2*k+D));



