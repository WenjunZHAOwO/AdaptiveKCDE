clear all;
close all;
rs = 5;
rng(rs);

path(path,'utils')

example = 6;
N = 800;
switch example
    case 0
        z_marginal = 'normal';
        zstar = [0.2,0.8];
%     case 1
%         z_marginal = 'uniform';
        
    case 2
        z_marginal = 'beta';
        zstar = [0.2, 0.7];
    case 6
        z_marginal = 'uniform';
        zstar = [0.3, 0.5];
       
end
weight = 0;
viz = true;
% zstar = [0.3, 0.5];
% nums_to_do = [1,2,3,5,6];
% nums_to_do = [1:3,5,6,8];
nums_to_do = 1:6;%[1,2,3,4,5,6];
num_eigs = 2;

center_x = true;

correction = false;

[x,z, fstar_true, z_type, range_star] = generate_data_2D(example, N, z_marginal, viz, zstar, rs);


if z_type == 'd'
    zstar = 1/2;
end





G = @(x,h,m) 1./(sqrt(2*pi*h.^2)).*exp( -1/2*( x - m').^2./h./h);



L = size(z,1);
hz = cell(1,6);
hx = cell(1,6);


if sum(1 == nums_to_do)
    % Proposal 1: benchmark: both rule-of-thumb
    [hx{1},hz{1},output1, lklhdval(1), KCDE_fun{1}, KCDE_fun_loo{1}] = find_both_h2(x,z, z_type,'none','constant',center_x); 
end

if sum(2 == nums_to_do)
    % Proposal 2: rule-of-thumb hz, variable hx
    [hx{2},hz{2},output2, lklhdval(2), KCDE_fun{2}, KCDE_fun_loo{2}] = find_both_h2(x,z, z_type,'none','variable',center_x);  
end

if sum(3 == nums_to_do)
    % Proposal 3: variable hx as a function of hz; hz constant chosen via LOO-CV
    [hx{3},hz{3},output3, lklhdval(3), KCDE_fun{3}, KCDE_fun_loo{3}] = find_both_h2(x,z, z_type,'z','variable',center_x); 
end


if sum(4 == nums_to_do)
    % Proposal 5: Solve for hz and then hx4 separately
    [hx{4},hz{4},output4, lklhdval(4), KCDE_fun{4}, KCDE_fun_loo{4}] = find_both_h2(x,z, z_type,'both','constant',center_x); 
end

if sum(5 == nums_to_do)
    % Proposal 6: Solve an optimization problem for two variables at once
    [hx{5}, hz{5},output5, lklhdval(5), KCDE_fun{5}, KCDE_fun_loo{5}] = find_both_h2(x,z, z_type,'both','variable',center_x); 
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

fprintf('finished solving bandwidths!')


%%


xl = [min( x(1,:))-1*std(x(1,:)),max( x(1,:))+1*std(x(1,:))];
yl = [min( x(2,:))-1*std(x(2,:)),max( x(2,:))+1*std(x(2,:))];

xstar1 = linspace(xl(1),xl(2),1e2);
xstar2 = linspace(yl(1),yl(2),1e2);
[Xstar1,Xstar2] = meshgrid(xstar1,xstar2);
fstar = cell(1,6);
score = zeros(1,6);
score_loo = zeros(1,6);

calc_score = false;

for i = 1:size(hz,2)
    i
    if size(hz{i},1)>0
        for j = 1:size(zstar,2)
%             fstar{i,j} = KCDE_fun{i}([Xstar1(:)'; Xstar2(:)'], zstar(:,j));
            
            if i == 6
                fstar{i,j} = KCDE_fun{6}{j}([Xstar1(:)'; Xstar2(:)'], zstar(:,j));
            else
                fstar{i,j} = KCDE_fun{i}([Xstar1(:)'; Xstar2(:)'], zstar(:,j));
            end
        end
        if calc_score
            if i ~= 6
                score(i) = 0;
                score_loo(i) = 0;
                score(i) = mean(diag( log(  KCDE_fun{i}(x, z))));

                for j = 1:N
                    %score(i) = score(i) + log(  KCDE_fun{i}(x(:,j), z(j)));
                    score_loo(i) = score_loo(i) + log(  KCDE_fun_loo{i}(x(:,j), z(j),  j));
                end
                %score(i) = score(i)/N;
                score_loo(i) = score_loo(i)/N;
             else

                 for j = 1:N
                        j/N
                        [~, ~,~, ~, KCDE_fun6, KCDE_fun_loo6] = find_both_h2_ultimate(x,z,z_type,'both','variable',center_x, z(j)); 


                        score(i) = score(i) + log(  KCDE_fun6(x(:,j), z(j)));
                        score_loo(i) = score_loo(i) + log(  KCDE_fun_loo6(x(:,j), z(j), j));


                 end
                 score(i) = score(i)/N;
                 score_loo(i) = score_loo(i)/N;
            end
        end
    end
    
end

score_true = sum( log(fstar_true(x,z)))/N;

% legs{1} = ['rule-of-thumb hz & hx  (',num2str(score(1),'%.2f'),')'];
% legs{2} = ['rule-of-thumb hz, variable hx  (',num2str(score(2),'%.2f'),')'];
% legs{3} = ['ML constant hz, variable hx  (',num2str(score(3),'%.2f'),')'];
% legs{4} = ['ML constant hz, constant hx, both opt (',num2str(score(4),'%.2f'),')'];
% legs{5} = ['ML constant hz, variable hx (\alpha), both opt (',num2str(score(5),'%.2f'),')'];
% legs{6} = ['ML constant hz, ultimate hx (\alpha etc), both opt  (',num2str(score(6),'%.2f'),')'];
    
%%
for j = 1:size(zstar,2)
    
    figure;
    counter = 1;
    
    subplot(1,7,1)
%     counter = counter + 1;
    contourf(Xstar1,Xstar2,reshape(fstar_true([Xstar1(:)'; Xstar2(:)'], zstar(j)),1e2,1e2),10);
    % title({'truth',['(',num2str(score_true,'%.2f'),')']})
    title(['Truth with $z=',num2str(zstar(j)),'$'],'interpreter','latex')
    xlim(xl)
    ylim(yl)
    cl = caxis;
%             caxis([0.0000    0.03])
%             caxis([0.0000    0.3])
    set(gca,'fontsize',16)
    axis square
    
    for i = 1:6
        if size(hz{i},1)>0
            subplot(1, 7,counter+1)
            counter = counter + 1;
            contourf(Xstar1,Xstar2,reshape(fstar{i,j},1e2,1e2));
            title(['Approach ',num2str(i)],'interpreter','latex')
            xlim(xl)
            ylim(yl)
            caxis(cl);
%             caxis([0.0000    0.03])
%             caxis([0.0000    0.3])
            set(gca,'fontsize',16)
            axis square
            
        end
    end
    set(gcf,'position',[-68 275 1509 258])


    % set(gcf,'position',[1478 920 911 883])
end% set(gcf,'position',[1478 520 511 883])



score
score_loo
lklhdval/size(z,2)
%%
zgrid = linspace(0,1,50);
if z_type == 'd'
    zgrid = unique(z);
end
Hdist = zeros( size(hz,2), length(zgrid));

for i = 1:size(hz,2)
    i
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
                Hdist(i, j) = pqdistances(f_est, f_true, 2, 'Hellinger', range_star(zgrid(j)));
            end
            
        end
end

%%
linestyles = {'-','--','-.','-','-.',':'};
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

