function [x,z, fstar_true,z_type, range_star] = generate_data_2D(example, N, z_marginal, viz, zstar, rs)

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
        z = [randn(1,round(N/2))+2, randn(1,round(N/2))-2];
        
    case 'categorical'
        z = [zeros(1,round(N/3)), 1/2*ones(1,round(N/4)), ones(1,round(N - N/3 - N/4)), ];
        z_type = 'd';
        
    case 'beta'
        z = betarnd(2,3,1,N);
%         z = betarnd(0.5,0.5,1,N);
%         z = betarnd(1.2,1.2,1,N);
        
end
z = (z-min(z))/(max(z)-min(z)); % normalize z to [0,1]

switch example
    
    case 0
        f=@(z) 5*[-sin(2*pi*z);2 + (1-z).^2];
        g=@(z)  2*[ 0.3*cos(2*pi*z+pi/4) + 0.4; 0.3*cos(2*pi*z+pi*3/4) + 0.4  ];
        
        x= f(z) + randn(2,N).*g(z);  
        
%         xstar = f(zstar) + randn(2,1e4).*g(zstar);   
        fstar_true = @(x,zstar) normpdf(x(1,:),[1,0]*f(zstar),[1,0]*g(zstar)).*normpdf(x(2,:),[0,1]*f(zstar),[0,1]*g(zstar));
        range_star = @(zstar) [f(zstar) - 6*g(zstar), f(zstar) + 6*g(zstar)];
        
    case 1
          
         
        f=@(z) 5*[-sin(2*pi*z);2 + (1-z).^2];
        g=@(z)  10*[0.1*z + 0.03;0.1*z + 0.03];
        
        x= f(z) + randn(2,N).*g(z);  
        
%         xstar = f(zstar) + randn(2,1e4).*g(zstar);   
        fstar_true = @(x,zstar) normpdf(x(1,:),[1,0]*f(zstar),[1,0]*g(zstar)).*normpdf(x(2,:),[0,1]*f(zstar),[0,1]*g(zstar));
        range_star = @(zstar) [f(zstar) - 6*g(zstar), f(zstar) + 6*g(zstar)];
        
    case 2
        
             
          
        f=@(z) [-sin(2*pi*z); 4*(z-1/2) ];
        g=@(z)  0.2;
        
        x= f(z) + randn(2,N).*g(z);
        x = x.^3;
        
        fstar_true = @(x,zstar) normpdf(sign(x(1,:)).*power(abs(x(1,:)),1/3),[1,0]*f(zstar),g(zstar)).*...
            normpdf(sign(x(2,:)).*power(abs(x(2,:)),1/3),[0,1]*f(zstar),g(zstar)).*...
            1/9.*power(abs(x(1,:))+1E-4,-2/3).*power(abs(x(2,:))+1E-4,-2/3);
        
        range_star = @(zstar) [f(zstar) - 6*g(zstar), f(zstar) + 6*g(zstar)].^3;
        
        
    case 3
       
        x = gamrnd(2,1,[2,N]);
%         xstar = gamrnd(2,1,[2,1e4]);
        
        
        fstar_true = @(x,zstar) gampdf(x(1,:),2,1).*gampdf(x(2,:),2,1);
        
    case 4
       
        x = cos(4*[z;1-z]) + 1/2*[.1;0.1].*randn(2,N);
        
        xstar = cos(4*[zstar;1-zstar]) + 1/2*[.1;0.1].*randn(2,1e4); 
        
        fstar_true = @(x,zstar) normpdf(x(1,:),cos(4*zstar),1/2*0.1).*normpdf(x(2,:),cos(4*(1-zstar)),1/2*0.1);

    case 5
        
        x = 2*[cos(2*pi*z);sin(2*pi*z)] + randn(2,N);
        idx1 = randperm(N); idx2 = randperm(N);
        idx1 = idx1(1:round(N/2)); idx2 = idx2(1:round(N/2));
        
        x(1,idx1) = -x(1,idx1);  x(2,idx2) = -x(2,idx2);
        
        xstar = 2*[cos(2*pi*zstar);sin(2*pi*zstar)] + randn(2,N);
        xstar(1,idx1) = -xstar(1,idx1); xstar(2,idx2) = -xstar(2,idx2);
        
        fstar_true = @(x,zstar) (1/2*normpdf(x(1,:),2*cos(2*pi*zstar),1)+1/2*normpdf(x(1,:),-2*cos(2*pi*zstar),1)).*...
             (1/2*normpdf(x(2,:),2*sin(2*pi*zstar),1)+1/2*normpdf(x(2,:),-2*sin(2*pi*zstar),1));

    case 6
        w = pi/4;
        f = @(z) 4*(z-z.^2).*[cos(2*pi*z+w);sin(2*pi*z+w)];
        g = @(z) [(z+1/16)/2;(1-z+1/16)/2];
        x = f(z) + g(z).*randn(2,N);
        idx1 = randperm(N); idx2 = randperm(N);
        idx1 = idx1(1:round(N/2)); idx2 = idx2(1:round(N/2));
        
        x(1,idx1) = -x(1,idx1)  ;  x(2,idx2) = -x(2,idx2) ;
        
%         xstar = cell(1,size(zstar,2));
%         for i = 1:size(zstar,2)
%             xstar{i} = 2*[cos(2*pi*zstar(i));sin(2*pi*zstar(i))] + 2*[(zstar(i)+1/16)/2;(zstar(i)+1/16)/2].*randn(2,N);
%             xstar{i}(1,idx1) = -xstar{i}(1,idx1) - zstar; xstar{i}(2,idx2) = -xstar{i}(2,idx2) - zstar(i);
%         end
        
        fstar_true = @(x,zstar) (1/2*normpdf(x(1,:),(-zstar.^2+zstar).*4.*cos(2*pi*zstar+w),(zstar+1/16)/2)+1/2*normpdf(x(1,:),(-zstar.^2+zstar).*(-4.*cos(2*pi*zstar+w)),(zstar+1/16)/2)).*...
             (1/2*normpdf(x(2,:),(-zstar.^2+zstar).*4.*sin(2*pi*zstar+w),(1-zstar+1/16)/2)+1/2*normpdf(x(2,:),(-zstar.^2+zstar).*(-4*sin(2*pi*zstar+w)),(1-zstar+1/16)/2));

        range_star = @(zstar) [-abs(f(zstar)) - 6*abs(g(zstar)), abs(f(zstar)) + 6*abs(g(zstar))];
        
end

if viz

    figure;
    idx = knnsearch(z',zstar');


    scatter(x(1,:),x(2,:),[],z,'filled')
    hold on,plot(x(1,idx),x(2,idx),'*','markersize',10,'linewidth',2)
    xlabel('$x_1$','interpreter','latex')
    ylabel('$x_2$','interpreter','latex')
    title('$x$ colored by $z$','interpreter','latex')

    set(gca,'fontsize',24)
    colorbar()
end
