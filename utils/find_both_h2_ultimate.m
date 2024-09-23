function [hx,hz,output,lklhdval, KCDE_fun, KCDE_fun_loo] = find_both_h2_ultimate(x, z, z_type, type,  hx_type, center_x, zstar, hz_weight, hz5)


N = size(x,2);
D = size(x,1);
output = [];
lklhdval = 0;

L = size(z,1);

if nargin < 4
    type = 'both';
end

if nargin < 5
    hx_type = 'variable';
end

if nargin < 6
    center_x = true;
end







options = optimset('MaxFunEvals',5e3);





% 
G = @(x,h,m) 1./(sqrt(2*pi*h.^2)).*exp( -1/2*( x - m').^2./h./h);


Gz = @(z,hz,mz) Gz_param(z, hz, mz, z_type);


optim_flag = 0;
% switch D
%     case 1
%         def_KCDE_funs;
%     case 2
%         def_KCDE_funs;%_2d;
% end
def_KCDE_funs;
% def_KCDE_funs_highDz;
hxz_param;%_highD_z;


if nargin < 8 || sum(hz_weight) == 0
    hz_weight = beta_init;
end


w = Gz(zstar, hz_weight, z)';



%loglklhd = @(hx, hz) scores_KCDE2(x,z,@(xstar, zstar, idx) KCDE_fun_loo(xstar, zstar, hx, hz, idx), true );
loglklhd = @(hx, hz) scores_KCDE2(x,z,[], true, true, KCDE_fun_loo_fast(hx, hz), w );
  


switch type
    
    case 'none'
        
        % hz = @(zstar) hz_fun(beta_init,zstar);
        hz = beta_init;
        switch hx_type
            case 'constant'
                hx = @(zstar) hx_fun( gamma_init,  zstar, hz, hz );
                

            case 'variable'
                hx = @(zstar) hx_fun( 1, zstar, hz, hz );
        end
        
        
    case 'z'
        hxfun_new = @(beta, zstar) hx_fun(1, zstar, hz_fun(beta,z) , beta );
        obj = @(beta) -loglklhd( @(z) hxfun_new(beta,z), beta );
        
        [output, outputval] = fminsearch(obj, beta_init ,options );
        output = output';
        %hz = @(zstar) hz_fun(output', zstar);
        hz = output;
        hx = @(zstar) hxfun_new(output', zstar);
            
        
        
    case 'both'

        hxfun_new = @(gamma, beta, zstar) hx_fun(gamma, zstar, hz_fun(beta,z) , beta );
        obj = @(gamma, beta) -loglklhd( @(z) hxfun_new(gamma, beta,z), beta );
 
       
        if nargin < 9 %|| hz5 == 0
            
            [output, outputval] = fminsearch(@(params) obj(params(1:gamma_D), params(gamma_D+1:end)'), [gamma_init, beta_init'] ,options );


            beta = output(gamma_D+1:end)';
            gamma = output(1:gamma_D);
        else

            [output, outputval] = fminsearch(@(params) obj( abs(params(1:gamma_D)), hz5), [gamma_init] ,options );
        
        
                beta = hz5;
                gamma = output(1:gamma_D);
                output = [output, hz5];
        end
        
        
        
        hz = beta;
        hx = @(zstar) hxfun_new(gamma, beta, zstar);
       
        
      
        
end




optim_flag = 1;
% switch D
%     case 1
%         def_KCDE_funs;
%     case 2
%         def_KCDE_funs_2d;
% end
def_KCDE_funs;%_highDz;

lklhdval = scores_KCDE2(x,z,KCDE_fun_loo, true )*N;
