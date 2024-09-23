% Main 1D
clear all;
close all;

path(path,'utils')
% Random seeds used in paper: 6,5,5,5 for example 0, 2, 6, 11
rs = 5;
% rng(rs);
example = 0;

N = 1000;
N_rs = 1;
weight = 0.0;
% N = 1600;
switch example
    case 0
        z_marginal = 'uniform';
        zstar = [0.21,0.75];
        rs = 6;
        xl = [-0.8,0.8];
        yl = [0,6];
    case 2
        z_marginal = 'beta';
        zstar = [0.4,0.8];
        rs = 5;
        xl = [-2.5,3.5];
        yl = [0,2.8];
    case 6
        z_marginal = 'uniform';
        zstar = [0.2,0.9];
        rs = 5;
        xl = [-1.3,3.8];
        yl = [0,3.5];
    case 10
        z_marginal = 'categorical2';
        zstar = [0.0,0.75];
        rs = 5;
        xl = [-5,10];
        yl = [0,0.3];
    case 12
        1
        z_marginal = 'uniform';
        zstar = [0.0,0.75];
        rs = 5;
        xl = [-22,22];
        yl = [0,1];
    case 13
        z_marginal = 'uniform';
        zstar = [0.49,0.51];
        rs = 5;
        xl = [-22,22];
        yl = [0,1];
    case 14
        z_marginal = 'uniform';
        zstar = [0.2,0.75];
        rs = 5;
        xl = [-10,5];
        yl = [0,2];
end


a = 2.5;
viz = true;


% zstar = [0, 0.75];

nums_to_do = [1:6];
% nums_to_do = 1:5;
% nums_to_do = [1,2,3,5,6,9];

center_x = true;


[x,z, xstar, fstar_true, z_type, range_star] = generate_data_1D(example, N, z_marginal, viz, zstar, rs, a);
x = x(1:N);
z = z(1:N);
if z_type == 'd'
    zstar = unique(z);
end


G = @(x,h,m) 1./(sqrt(2*pi*h.^2)).*exp( -1/2*( x - m').^2./h./h);

L = size(z,1);
hz = cell(1,6);
hx = cell(1,6);
lklhdval = zeros(1,6);


if sum(1 == nums_to_do)
    % Proposal 1: benchmark: both rule-of-thumb
    [hx{1},hz{1},output1, lklhdval(1), KCDE_fun{1}, KCDE_fun_loo{1}] = find_both_h2(x,z,z_type,'none','constant',center_x); 
end

if sum(2 == nums_to_do)
    % Proposal 2: rule-of-thumb hz, variable hx
    [hx{2},hz{2},output2, lklhdval(2), KCDE_fun{2}, KCDE_fun_loo{2}] = find_both_h2(x,z,z_type,'none','variable',center_x); 
end

if sum(3 == nums_to_do)
    % Proposal 3: variable hx as a function of hz; hz constant chosen via LOO-CV
    [hx{3},hz{3},output3, lklhdval(3), KCDE_fun{3}, KCDE_fun_loo{3}] = find_both_h2(x,z,z_type,'z','variable',center_x); 
end



if sum(4 == nums_to_do)
    % Proposal 4: Solve for hz and then hx4 separately
    [hx{4},hz{4},output4, lklhdval(4), KCDE_fun{4}, KCDE_fun_loo{4}] = find_both_h2(x,z,z_type,'both','constant',center_x); 
end

if sum(5 == nums_to_do)
    % Proposal 5: Solve an optimization problem for two variables at once
    [hx{5}, hz{5},output5, lklhdval(5), KCDE_fun{5}, KCDE_fun_loo{5}] = find_both_h2(x,z,z_type,'both','variable',center_x); 
end


if sum(6 == nums_to_do)
    % Proposal 6: Solve an optimization problem for each zstar
    hx{6} = cell(1,length(zstar));
    hz{6} = cell(1,length(zstar));
    KCDE_fun{6} = cell(1,length(zstar));
    KCDE_fun_loo{6} = cell(1,length(zstar));
    for i = 1:length(zstar)
        [hx{6}{i}, hz{6}{i},output6, lklhdval(6,i), KCDE_fun{6}{i}, KCDE_fun_loo{6}{i}] = find_both_h2_ultimate(x,z,z_type,'both','variable',center_x, zstar(i)); 

    end
end






score_true = sum( log(fstar_true(x,z)))/N;
linestyles = {'-','--','-.','-','-.',':'};


%%
calc_score = false;
if calc_score
    for i = 1:size(hz,2)
            if size(hz{i},1)>0
                score(i) = 0;
                score_loo(i) = 0;
                if i ~= 6
                    for j = 1:N
                        score(i) = score(i) + log(  KCDE_fun{i}(x(j), z(j)));
                        score_loo(i) = score_loo(i) + log(  KCDE_fun_loo{i}(x(j), z(j), j));

                    end
                else
                    for j = 1:N
                        j/N
                        [~, ~,~, ~, KCDE_fun6, KCDE_fun_loo6] = find_both_h2_ultimate(x,z,z_type,'both','variable',center_x, z(j)); 


                        score(i) = score(i) + log(  KCDE_fun6(x(j), z(j)));
                        score_loo(i) = score_loo(i) + log(  KCDE_fun_loo6(x(j), z(j), j));


                    end
            %                 
            %                 
                end




                score(i) = score(i)/N;
                score_loo(i) = score_loo(i)/N;    
            end
    end
end

%%
% figure;
for l = 1:size(zstar,2)
    
    %xl = [min(x)-std(x), max(x)+std(x)];
    xgrid = linspace(xl(1),xl(2),1e3);
    zstarl = zstar(l);%*ones(1,length(xgrid));

    fstar = cell(1,size(hz,2)); 
    

    figure;


    for i = 1:size(hz,2)
        if size(hz{i},1)>0
            if i == 6
                fstar{i} = KCDE_fun{6}{l}(xgrid, zstarl);
            else
                fstar{i} = KCDE_fun{i}(xgrid, zstarl);
            end
            
                
                
            
            
            
            plot(xgrid,fstar{i},linestyles{i},'linewidth',2), hold on,
            
        end
    end
    
    

    plot(xgrid,fstar_true(xgrid,zstar(l)),'k','linewidth',2)
    
    set(gca,'fontsize',16)
    ylim(yl);
    xlim(xl);
%     legs{1} = ['rule-of-thumb hz & hx  (',num2str(score(1),'%.2f'),')'];
%     legs{2} = ['rule-of-thumb hz, variable hx  (',num2str(score(2),'%.2f'),')'];
%     legs{3} = ['ML constant hz, variable hx  (',num2str(score(3),'%.2f'),')'];
%     legs{4} = ['ML constant hz, constant hx, both optimized (',num2str(score(4),'%.2f'),')'];
%     legs{5} = ['ML constant hz, variable hx (\alpha), both optimized  (',num2str(score(5),'%.2f'),')'];
%     legs{6} = ['ML constant hz, ultimate hx (\alpha etc), both optimized  (',num2str(score(6),'%.2f'),')'];
%     
%     legs_used = cell(1,length(nums_to_do));
%     counter = 1;
%     for i = 1:6
%         if sum(i == nums_to_do)
%             legs_used{counter} = legs{i};
%             counter = counter + 1;
%         end
%     end
%     leg = legend(legs_used);
title(['$\rho(x|$',num2str(zstar(l)),'$)$'],'interpreter','latex')
legend('Approach 1','Approach 2','Approach 3','Approach 4','Approach 5','Approach 6')
% set(leg,'location','southoutside','fontsize',14)
xlabel('$x$','interpreter','latex')





% set(gcf,'position',[-647 1047 611 618]);

set(gcf,'position',[270 1230 632 328]);
end

if calc_score


    score
    score_loo
    score_true
    lklhdval/N
end

%%

zgrid = sort(z);
switch z_type
    case 'c'
        zgrid = linspace(0,1,50);
    case 'd'
        zgrid = unique(z);
end
Hdist = zeros( size(hz,2), length(zgrid));
ISEdist = zeros( size(hz,2), length(zgrid));
for i = 1:size(hz,2)
        if size(hz{i},1)>0
            for j = 1:length(zgrid)
                if i == 6
                    j/length(zgrid)
                    [~, ~, output6, ~, KCDE_fun6, KCDE_fun_loo6] = find_both_h2_ultimate(x,z,z_type,'both','variable',center_x, zgrid(j)); 
                    hz6(j) = abs(output6(2));
                    alpha6(j) = abs(output6(1));
                    f_est = @(xi) KCDE_fun6(xi, zgrid(j));
                    
                    
                    figure(19);
                    subplot(1,2,1)
                    switch z_type
                        case 'c'
                            scatter(zgrid(j), hz6(j),'k','filled'); hold on,
                            hold on,plot(linspace(0,1,100),output5(2)*ones(1,100),'b-.','linewidth',2)
                        case 'd'
                            scatter(zgrid(j),  (1/2+1/2*cos(2*pi*hz6(j)).^2),'k','filled'); hold on,
                            hold on,plot(linspace(0,1,100),(1/2+1/2*cos(2*pi*output5(2)).^2)*ones(1,100),'b-.','linewidth',2)
                            
                    end
                    
                    set(gca,'fontsize',16)
                    xlabel('$z$','interpreter','latex')
                    title('$h_z(z)$','interpreter','latex')
                    xlim([min(z),max(z)]);
                    
                    

                    subplot(1,2,2)
                    scatter(zgrid(j), alpha6(j),'k','filled'); hold on,
                    
                    set(gca,'fontsize',16)
                    xlabel('$z$','interpreter','latex')
                    title('$\alpha(z)$','interpreter','latex')
                    xlim([min(z),max(z)]);
                    set(gcf,'position',[21 452 853 286])
                    hold on,plot(linspace(0,1,100),output5(1)*ones(1,100),'b-.','linewidth',2)

                    
                else
                    f_est = @(xi) KCDE_fun{i}(xi, zgrid(j));
                end
                f_true = @(xi) fstar_true(xi,zgrid(j));
                Hdist(i, j) = pqdistances(f_est, f_true, 1,'Hellinger',  range_star(zgrid(j)) );
%                 ISEdist(i, j) = pqdistances(f_est, f_true, 1,'ISE',range_star(zgrid(j)));
                
            end
            
        end
end

mean(Hdist')

%%
figure;
for i = 1:size(hz,2)
    if size(hz{i},1)>0
        
        plot( zgrid, Hdist(i,:),linestyles{i},'linewidth',2), hold on
    end
end
legs{1} = ['rule-of-thumb hz & hx  (',num2str(mean(Hdist(1,:)),'%.2f'),')'];
legs{2} = ['rule-of-thumb hz, variable hx  (',num2str(mean(Hdist(2,:)),'%.2f'),')'];
legs{3} = ['ML constant hz, variable hx  (',num2str(mean(Hdist(3,:)),'%.2f'),')'];
legs{4} = ['ML constant hz, constant hx, both optimized (',num2str(mean(Hdist(4,:)),'%.2f'),')'];
legs{5} = ['ML constant hz, variable hx (\alpha), both optimized  (',num2str(mean(Hdist(5,:)),'%.2f'),')'];
legs{6} = ['ML constant hz, ultimate hx (\alpha etc), both optimized  (',num2str(mean(Hdist(6,:)),'%.2f'),')'];
    

legs_used = cell(1,length(nums_to_do));
counter = 1;
for i = 1:6
    if sum(i == nums_to_do)
        legs_used{counter} = legs{i};
        counter = counter + 1;
    end
end
% leg = legend(legs_used);
% 
% set(leg,'location','southoutside','fontsize',14)
xlabel('$z$','interpreter','latex')
title('squared Hellinger distance vs $z$','interpreter','latex')
set(gca,'fontsize',16)


legend('Approach 1','Approach 2','Approach 3','Approach 4','Approach 5','Approach 6')
% set(leg,'location','southoutside','fontsize',14)
xlabel('$x$','interpreter','latex')


set(gcf,'position',[-647 1047 611 618]);

set(gcf,'position',[270 1230 632 328]);

%%
% if strcmp(z_marginal, 'categorical')
%     
%     
%     
%     hx_vals = zeros( size(hz,2), length(zgrid));
% 
%     for i = 1:size(hz,2)
%             if size(hz{i},1)>0
%                 for j = 1:length(zgrid)
%                     hx_vals(i,j) = hx{i}(zgrid(j));
%                 end
% 
%             end
%     end
%     
%     
%     
% % end
% 
% figure;
% for i = 1:size(hz,2)
%     if size(hz{i},1)>0
%         
%         plot( zgrid, hx_vals(i,:),'linewidth',2), hold on
%     end
% end
% legs{1} = ['rule-of-thumb hz & hx  (',num2str(mean(hx_vals(1,:)),'%.2f'),')'];
% legs{2} = ['rule-of-thumb hz, variable hx  (',num2str(mean(hx_vals(2,:)),'%.2f'),')'];
% legs{3} = ['ML constant hz, variable hx  (',num2str(mean(hx_vals(3,:)),'%.2f'),')'];
% legs{4} = ['ML constant hz, constant hx, both optimized (',num2str(mean(hx_vals(4,:)),'%.2f'),')'];
% legs{5} = ['ML constant hz, variable hx (\alpha), both optimized  (',num2str(mean(hx_vals(5,:)),'%.2f'),')'];
% legs{6} = ['ML constant hz, ultimate hx (\alpha etc), both optimized  (',num2str(mean(hx_vals(6,:)),'%.2f'),')'];
%   
% legs_used = cell(1,length(nums_to_do));
% counter = 1;
% for i = 1:6
%     if sum(i == nums_to_do)
%         legs_used{counter} = legs{i};
%         counter = counter + 1;
%     end
% end
% leg = legend(legs_used);
% 
% set(leg,'location','southoutside','fontsize',14)
% xlabel('$z$','interpreter','latex')
% title('$hx(z)$ vs $z$','interpreter','latex')
% set(gca,'fontsize',16)
% 
% mean(Hdist')
% 
% 
% %%
% 
xl = [min(x)-std(x), max(x)+std(x)];
xgrid = linspace(xl(1),xl(2),100);

zgrid = sort(z);
switch z_type
    case 'c'
        zgrid = linspace(0,1,20);
    case 'd'
        zgrid = unique(z);
end
% xgrid = linspace(min(xgrid),max(xgrid),100);

ftrue_array = zeros(length(zgrid),length(xgrid));
for i = 1:length(zgrid)
    
   ftrue_array(i,:) = fstar_true(xgrid, zgrid(i)); 
    
end

figure;
% if z_type == 'c'
%     [X,Z] = meshgrid(xgrid, zgrid);
%     s=surf(X,Z,ftrue_array,'FaceAlpha',0.4);
%     s.EdgeColor='none';
%     
% else
    for i = 1:length(zgrid)
        p=fill3( xgrid, zgrid(i)+0*xgrid, ftrue_array(i,:), zgrid(i) );hold on,
        p.FaceAlpha = 0.4;
        p.EdgeColor = 'none';
    end
% %     [X,Z] = meshgrid(xgrid, zgrid);
% %     s=surf(X,Z,ftrue_array,'FaceAlpha',0.4);
%     %s.EdgeColor='none';
% end
xlabel('$x$','interpreter','latex')
ylabel('$z$','interpreter','latex')
set(gca,'fontsize',16)
title('$\rho^{true}(x|z)$','interpreter','latex')


hold on,
scatter(x,z,10,'k','filled')
xl = xlim;
if z_type == 'csdfsdf'
    for l = 1:size(zstar,2)
        plot(linspace(xl(1),xl(2),N),zstar(l)*ones(1,N),'k-.','linewidth',2)
    end
end
xlim(xl);
box on
grid on

