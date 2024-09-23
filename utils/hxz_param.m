gamma_D = 1;
beta_init = zeros(L,1);
Lc = sum(z_type == 'c')+2*sum(z_type=='p');
for l = 1:L
    switch z_type(l)
        case 'c'
            
            beta_init(l) = (4/(Lc+2))^(1/(Lc+4))*std(z(l,:)).*N.^(-1/(4+Lc));
        case 'd'
            beta_init(l) = 1;
        case 'p'
            beta_init(l) = (4/(Lc+2))^(1/(Lc+4))*std([cos(2*pi*z(l,:)),sin(2*pi*z(l,:))]).*N.^(-1/(4+Lc));
    end
end

% if pseudo_cont && strcmp(type, 'none')
%    beta_init =  0.01*ones(L,1);
% end


switch hx_type
    case 'variable'
        
        
        hx_actual = @(zstar, hz, hzstar) rule_of_thumb_x(zstar, hz, x,z,Gz,hzstar);
        
        hx_fun = @(gamma, zstar, hz, hzstar) gamma*hx_actual(zstar, hz, hzstar);
        
            
        gamma_init = 1;
        
    case 'constant'
        hx_fun = @(gamma, zstar, hz, hzstar) gamma'.*ones(1,size(zstar,2));
        gamma_init = (4/(D+2))^(1/(D+4))*std(x').*N.^(-1/(4+D));
            
        gamma_D = D;
        
        
    case 'ultimate'
        
        hx_actual = @(zstar, hz, hzstar) rule_of_thumb_x(zstar, hz, x,z,Gz,hzstar);
        
        hx_fun = @(gamma, zstar, hz, hzstar) (gamma*[ones(1,size(zstar,2));zstar] ).*hx_actual(zstar, hz, hzstar);
        
            
        gamma_init = [1,zeros(1,size(z,1))];
        gamma_D = length(gamma_init);
    case 'oracle'
        hx_actual = @(zstar, hz, hzstar,k) rule_of_thumb_x(zstar, hz, x,z,Gz,hzstar,k);
        
        hx_fun = @(gamma, zstar, hz, hzstar, k) (gamma*[ones(1,size(zstar,2));zstar] ).*hx_actual(zstar, hz, hzstar,k);
        
            
        gamma_init = [1,zeros(1,size(z,1))];
        gamma_D = length(gamma_init);
        k_init = 2;
%     case 'faithful'
%         
%         hx_actual = @(zstar, hz, hzstar) rule_of_thumb_x(zstar, hz, x,z,Gz,hzstar);
%         
%         %hx_fun = @(gamma, zstar, hz, hzstar) (gamma(1)+gamma(2)*tanh(gamma(3:end)*[0*ones(1,size(zstar,2));(zstar-mean(z')')./std(z')']) ).*hx_actual(zstar, hz, hzstar);
%         hx_fun = @(gamma, zstar, hz, hzstar) (gamma(1)+gamma(2)^2*tanh(gamma(3)*(zstar-70)) ).*hx_actual(zstar, hz, hzstar);
%         
%         hx_fun = @(gamma, zstar, hz, hzstar) (gamma(1)*(zstar<=68) + gamma(2)*(zstar>68)).*hx_actual(zstar, hz, hzstar);
%         
%         gamma_init = [1,1];   
% %         gamma_init = [1, 1, 0 ];
% %         gamma_init = [-0.1671  -43.5850   69.7195    5.5132];
%         gamma_D = length(gamma_init);
        
%     case 'faithful2'
%         
%         hx_actual = @(zstar, hz, hzstar) rule_of_thumb_x(zstar, hz, x,z,Gz,hzstar);
%         
%         %hx_fun = @(gamma, zstar, hz, hzstar) (gamma(1)+gamma(2)*tanh(gamma(3:end)*[0*ones(1,size(zstar,2));(zstar-mean(z')')./std(z')']) ).*hx_actual(zstar, hz, hzstar);
%         hx_fun = @(gamma, zstar, hz, hzstar) (gamma(1)+gamma(2)^2*tanh(gamma(3)*(zstar-70)) ).*hx_actual(zstar, hz, hzstar);
%         
% %         hx_fun = @(gamma, zstar, hz, hzstar) (gamma(1)*(zstar<=68) + gamma(2)*(zstar>68)).*hx_actual(zstar, hz, hzstar);
%         
%         gamma_init = [1,1,0];   
% %         gamma_init = [1, 1, 0 ];
% %         gamma_init = [-0.1671  -43.5850   69.7195    5.5132];
%         gamma_D = length(gamma_init);
%         
        
    
               

end


