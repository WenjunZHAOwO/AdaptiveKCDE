function [x,z, xstar, fstar_true, z_type, range_star] = generate_data_1D(example, N, z_marginal, viz, zstar, rs, a)

if nargin < 6
    rs = 5;
end
if nargin < 7
    a = 0;
end
rng(rs)
z_type = 'c';
switch z_marginal % try different marginal distributions for z...
    case 'uniform'
        z = rand(1,N);
        
    case 'normal'
        z = randn(1,N);
       
    case 'skewed'
        z = rand(1,N).^3;
        
    case 'bimodal'
        z = [randn(1,round(N/3))+2, randn(1,round(N/3*2))-2];
        
    case 'beta'
        z = betarnd(2,3,1,N);
%         z = betarnd(0.5,0.5,1,N);
        
%         z = betarnd(0.5,1,1,N);
        
        
    case 'categorical'
        z = [zeros(1,round(N/20)), ones(1,round(N/10)), 2*ones(1,round(N/4)), 3*ones(1,round(N/5)), 4*ones(1,round(N/4) )];
        z = [z, 5*ones(1, round(N/10) )];
        z = [z, 6*ones(1, N-length(z)) ];
        z = 1-z;
        z_type = 'd';
        z = z(randperm(N));
        
   case 'categorical2'
        
        z = [zeros(1,round(N/10)), ones(1,round(6*N/20)), 2*ones(1,round(N/5)), 3*ones(1,round(N/20*6) )];
        z = [z, 4*ones(1, N-length(z)) ];
        z = 1-z;
        z_type = 'd';
        z = z(randperm(N));
        
end
z = (z-min(z')')/(max(z')-min(z'))'; % normalize z to [0,1]
% zstar = quantile(z,0.2);
% zstar = 1.2%0.99%[0.1,0.4, 0.9];
% zstar = 0.9;

xstar = zeros(1,N);

switch example
    
    case 0
        
        f=@(z) 0*z;
        g=@(z)  0.5*cos(4*pi*z + pi/4) + 0.55;
        
        x= f(z) + randn(1,N).*g(z);        
        for l = 1:size(zstar,2)      
            xstar(l,:) = f(zstar(l)) + randn(1,N).*g(zstar(l));   
        end
        fstar_true = @(x,zstar) normpdf(x,f(zstar),g(zstar));
        range_star = @(zstar) [f(zstar) - 6*g(zstar), f(zstar) + 6*g(zstar)];
        
    
    case 1
          
        
        
        f=@(z) -sin(2*pi*z);
        g=@(z)  0.5*z + 0.05;
        
        x= f(z) + randn(1,N).*g(z);        
        for l = 1:size(zstar,2)      
            xstar(l,:) = f(zstar(l)) + randn(1,N).*g(zstar(l));   
        end
        fstar_true = @(x,zstar) normpdf(x,f(zstar),g(zstar));
        range_star = @(zstar) [f(zstar) - 6*g(zstar), f(zstar) + 6*g(zstar)];
        
        
        
    case 2
        
        
        f=@(z) -sin(2*pi*z);
        g=@(z)  0.2 + 0*z;
        
        x= f(z) + randn(1,N).*g(z);
        x = x.^3;
              
        for l = 1:size(zstar,2)
            xstar(l,:) = f(zstar(l)) + randn(1,N).*g(zstar(l));   

            xstar(l,:) = xstar(l,:).^3;
        end
        
        fstar_true = @(x,zstar) normpdf( sign(x).*power(abs(x),1/3),f(zstar),g(zstar)).*(1/3.*power(abs(x) +1E-3,-2/3));
        range_star = @(zstar) [f(zstar) - 6*g(zstar), f(zstar) + 6*g(zstar)].^3;
        
        
    case 3
       
        x = trnd( 3*z + 1.1, 1, N); 
        
        xstar = [];
        fstar_true = @(x,zstar) tpdf( x, 3*zstar + 1.1);
        
        range_star = @(zstar) [-6,6];
        
        
     
    case 4
        x = 5*(tanh( 2*( z-1/2) )).*randn(1,N);
        
        
        for l = 1:size(zstar,2)
        
            xstar(l,:) = 5*(tanh( 2*( zstar(l)-1/2) )).*randn(1,N);
            
        end
        
        fstar_true = @(x,zstar)  normpdf(x,0, 5*abs((tanh( 2*( zstar-1/2) ))) );
    
    case 5
        
        a = rand(1,N);
        w = randn(1,N);
               
        x = (a<1/2).*(z<1/2).*(2+tanh(5*z-5/2)) + (a>1/2).*(-tanh(5*z-5/2)) +...
            + (a<1/2).*(z>1/2).*(tanh(5*z-5/2)-2);
        x = x + w*0.15;
        
               
        for l = 1:size(zstar,2)       
            xstar(l,:) = (a<1/2).*(zstar(l)<1/2).*(2+tanh(5*zstar(l)-5/2)) + (a>1/2).*(-tanh(5*zstar(l)-5/2)) +...
                + (a<1/2).*(zstar(l)>1/2).*(tanh(5*zstar(l)-5/2)-2);
            xstar(l,:) = xstar(l,:)+w*0.15;
        end
        
        
        fstar_true = @(x,zstar)  1/2*normpdf(x,-tanh(5*zstar-5/2),0.15) + ...
            + 1/2*normpdf( x, (zstar<1/2).*(2+tanh(5*zstar-5/2)) + (zstar>1/2).*(tanh(5*zstar-5/2)-2),0.15);
        
    case 6
        
        a = rand(1,N);
        w = .2*randn(1,N);
        
        g = @(z) z+1/4;
               
        x = (a <  0.3 ).*w + (a > 0.3 ).*(2 + w ) ;
        x = x.*g(z);
        
        
               
        xstar = [];
        
        
        
        fstar_true = @(x,zstar)  0.3.*normpdf(x,0,.2*g(zstar))+...
            + 0.7.*normpdf(x,2*g(zstar),.2*g(zstar));
        
        range_star = @(zstar) [ -6*.2*g(zstar), 2*g(zstar)+6*.2*g(zstar) ];
        
        
    case 7
        
        
               
        x = (0.5*abs(cos(2*pi*z))+0.05).*randn(1,N);
               
        for l = 1:size(zstar,2)       
            xstar(l,:) = (0.5*abs(cos(2*pi*zstar(l))) + 0.05).*randn(1,N);
        end
        
         
        fstar_true = @(x,zstar)  normpdf( x,0, abs(0.5*cos(2*pi*zstar))+0.05  );
          
    case 8
        
        x = randn(1,N);
        
        ra = rand(1,N);
        x = x + (-0.4*x + a*(1-z)).*(ra<1/8) + (-0.4*x -a*(1-z)).*(ra>7/8);
        
        for l = 1:size(zstar,2)
        
            xstar(l,:) = randn(1,N);
            
            xstar(l,:) = xstar(l,:) + (-0.4*xstar(l,:) + a*(1-zstar(l))).*(ra<1/8) +(-0.4*xstar(l,:) -a*(1-zstar(l))).*(ra>7/8);
            
        end
        
        
        
        fstar_true = @(x,zstar) 1/8.*normpdf(x,-a*(1-zstar),0.6) + 1/8.*normpdf(x,a*(1-zstar),0.6) + 3/4*normpdf(x,0,1);
        
        
        
    case 9
        
        x = z + cos(pi*z) + 0.01*randn(1,N);
        for l = 1:size(zstar,2)       
            xstar(l,:) = zstar(l) + cos(pi*zstar(l)) + 0.01*randn(1,N);
        end
        fstar_true = @(x,zstar) normpdf(x,zstar+cos(pi*zstar),0.01);
        
    case 10
        
 	% designed for categorical, bimodal ver.
        
        x = randn(1,N);
        
        ra = rand(1,N);
        x = x + (-0.4.*x + a*(1-z)).*(ra<1/8) + (-0.4.*x -a*(1-z)).*(ra>7/8);
        
        g = @(z) ( z + 2)/3;
        x = x.*g(z);
        
        idx = randperm(N);
        x(idx(1:round(N/2))) = x(idx(1:round(N/2))) + 5;       
        for l = 1:size(zstar,2)
        
            xstar(l,:) = randn(1,N);
            
            xstar(l,:) = xstar(l,:) + (-0.4*xstar(l,:) + a*(1-zstar(l))).*(ra<1/8) +(-0.4*xstar(l,:) -a*(1-zstar(l))).*(ra>7/8);
            
        end
        fstar_true_first = @(x,zstar) 1/8.*normpdf(x,-a*(1-zstar).*g(zstar),0.6*g(zstar) ) + 1/8.*normpdf(x,a*(1-zstar).*g(zstar),0.6.*g(zstar)) +...
            + 3/4*normpdf(x,0,1.*g(zstar));
        fstar_true = @(x,zstar) 1/2*fstar_true_first(x,zstar) + 1/2*fstar_true_first(x-5,zstar);
        
        
        range_star = @(zstar) g(zstar).*[-a*(1-zstar)-1.8, a*(1-zstar) + 1.8];
        range_star = @(zstar) [-5,10];      
    case 11
        % designed for categorical
        
        x = randn(1,N);
        
        ra = rand(1,N);
        x = x + (-0.4.*x + a*(1-z)).*(ra<1/8) + (-0.4.*x -a*(1-z)).*(ra>7/8);
        
        g = @(z) ( z + 2)/3;
        x = x.*g(z);
        
        for l = 1:size(zstar,2)
        
            xstar(l,:) = randn(1,N);
            
            xstar(l,:) = xstar(l,:) + (-0.4*xstar(l,:) + a*(1-zstar(l))).*(ra<1/8) +(-0.4*xstar(l,:) -a*(1-zstar(l))).*(ra>7/8);
            
        end
        fstar_true = @(x,zstar) 1/8.*normpdf(x,-a*(1-zstar).*g(zstar),0.6*g(zstar) ) + 1/8.*normpdf(x,a*(1-zstar).*g(zstar),0.6.*g(zstar)) +...
            + 3/4*normpdf(x,0,1.*g(zstar));
        
        
        range_star = @(zstar) g(zstar).*[-a*(1-zstar)-1.8, a*(1-zstar) + 1.8];
        range_star = @(zstar) [-3.5,3.5];
         
    case 12
        
        % smoothness 1
        g=@(z) tanh(10*(z-0.5))+1.5;
        x = randn(1,N).*g(z) + 10*tanh(10*(z-0.5));
        fstar_true = @(x,zstar) normpdf( x, 10*tanh(10*(zstar-0.5)), g(zstar)) ;
        range_star = @(zstar) [-22,22];
        range_star = @(zstar) 10*tanh(10*(zstar-0.5)) + 6*[-g(zstar), g(zstar)];
         
        
    case 13
        
        % smoothness 2
        g = @(z) tanh(10*(z-0.5))+1.5;
        x = g(z).*randn(1,N) + 20*(heaviside(z-0.5)-0.5);
        fstar_true = @(x,zstar) normpdf( x, 20*(heaviside(zstar-0.5)-0.5), g(zstar) );
        range_star = @(zstar) 20*(heaviside(zstar-0.5)-0.5) + 6*[-g(zstar), g(zstar)];

	case 14

		g = @(z) 3*tanh(10*(z-0.5));%+1.5;
        %xs = rand(1,N);

        %x = sqrt(xs)-1 ; idxx = randperm(N); x(idxx(1:round(N/2))) = -x(idxx(1:round(N/2)));
        %x = x + g(z);g = @(z) tanh(10*(z-0.5))+1.5;
        g_var = @(z) z + 1/4;
        xs = rand(1,N);
        
        x = sqrt(xs)-1 ; idxx = randperm(N); x(idxx(1:round(N/2))) = -x(idxx(1:round(N/2)));
        x = g_var(z).*x + g(z);
        f = @(x) (1-abs(x)).*( abs(x)<1);
        fstar_true = @(x,zstar) f( (x - g(zstar))./g_var(zstar) )./g_var(zstar);
        range_star = @(zstar) [g(zstar) - 2 - g_var(zstar), g(zstar)+2+g_var(zstar)];

        %f = @(x) (1-abs(x)).*( abs(x)<1);
        %fstar_true = @(x,zstar) f( x - g(zstar) );
        %range_star = @(zstar) [g(zstar) - 2, g(zstar)+2];
end  


if viz
    
    figure;
    scatter(z,x)
    xlabel('$z$','interpreter','latex')
    ylabel('$x$','interpreter','latex')
    title('$x$ vs $z$','interpreter','latex')
    set(gca,'fontsize',24)
    yl = ylim;
    hold on,
    for l = 1:size(zstar,2)
        plot(zstar(l)*ones(1,N),linspace(yl(1),yl(2),N),'k-.','linewidth',2)
    end
    leg=legend('data','$z^*$');
    set(leg,'interpreter','latex')
    ylim(yl);
    set(gcf,'position',[-33 1160 560 420]);



end
