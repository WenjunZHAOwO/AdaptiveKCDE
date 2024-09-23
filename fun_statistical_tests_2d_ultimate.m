% Main 1D
function fun_statistical_tests_2d_ultimate(example, rs)


N = 800;
% N = 1600;
switch example
    case 0
        z_marginal = 'uniform';
        fileName = ['Gaussian2_',num2str(rs),'_out.txt'];
    case 1
        z_marginal = 'uniform';
        fileName = ['Gaussian_',num2str(rs),'_out.txt'];
    case 2
        z_marginal = 'beta';
        fileName = ['Skew_',num2str(rs),'_out.txt'];
    case 6
        z_marginal = 'uniform';
        fileName = ['Mtmd_',num2str(rs),'_out.txt'];
end

% z_marginal = 'beta';
% z_marginal = 'bimodal';
a = 0;
viz = false;
zstar = [0.45, 0.8];
zstar = [0.35,0.9];
zstar = [0,1/4,1/2,3/4,1];
zstar = [0,1];
zstar = [0];
nums_to_do = [1,2,3,4,5,6];

weight = 0

fileID = fopen(fileName,'a');

center_x = true;

N_rs = 1;
hdists_rs = zeros(6,N_rs);
score_rs = zeros(6,N_rs);
score_loo_rs = zeros(6,N_rs);
alphas_rs = zeros(1,N_rs);
score_cv_rs = zeros(6,N_rs);


    rs
    tic
    [x,z,  fstar_true, z_type, range_star] = generate_data_2D(example, N+200, z_marginal, viz, zstar, rs);
    x_cv = x(:, N+1:end);
    z_cv = z(:, N+1:end);
    x = x(:, 1:N);
    z = z(:, 1:N);
    
%     x_cv = min( x_cv, max(x')');
%     x_cv = max( x_cv, min(x')');

    G = @(x,h,m) 1./(sqrt(2*pi*h.^2)).*exp( -1/2*( x - m').^2./h./h);

    L = size(z,1);
    hz = cell(1,5);
    hx = cell(1,5);
    lklhdval = zeros(1,5);


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
        % Proposal 6: Solve an optimization problem for two variables at once
        [hx{5}, hz{5},output5, lklhdval(5), KCDE_fun{5}, KCDE_fun_loo{5}] = find_both_h2(x,z,z_type,'both','variable',center_x); 
        alphas_rs(rs) = output5(1);
    end





    %%
    % figure;
%     for l = 1:size(zstar,2)

%         xl = [min(x)-std(x), max(x)+std(x)];
%         xgrid = linspace(xl(1),xl(2),1e3);
%         zstarl = zstar(l);%*ones(1,length(xgrid));
% 
%         fstar = cell(1,size(hz,2)); 
%         score = zeros(1,size(hz,2));
%         score_loo = zeros(1,size(hz,2));
% 
%         figure;


        for i = 1:6
            if ismember(i,nums_to_do)
%                 fstar{i} = KCDE_fun{i}(xgrid, zstarl);
                score(i) = 0;
                score_loo(i) = 0;
                score_cv(i) = 0;
                
                if i ~= 6
                    score(i) = mean(diag( log(  KCDE_fun{i}(x, z) )));

                    for j = 1:N
                        %score(i) = score(i) + log(  KCDE_fun{i}(x(j), z(j)));
                        score_loo_temp(j) = log(  KCDE_fun_loo{i}(x(:,j), z(j), j));
                        %score_loo(i) = score_loo(i) + log(  KCDE_fun_loo{i}(x(j), z(j), j));

                    end

                    score_loo(i) = mean(score_loo_temp);  
                
                    score_cv(i) = mean(diag( log(  KCDE_fun{i}(x_cv, z_cv) )));
                    if score_cv(i) == -Inf
                        stop
                    end
                    
                else
                    
                    score_loo_temp = zeros(1,N);
                    score_temp = zeros(1,N);
                    for j = 1:N
                        fprintf(fileID,['Approach 6 Run ',num2str(rs), ', in-sample: ',num2str(j/N),'\n'])
                        [~, ~,~, ~, KCDE_fun6, KCDE_fun_loo6] = find_both_h2_ultimate(x,z,z_type,'both','variable',center_x, z(j),weight); %                         score(i) = score(i) + log(  KCDE_fun6(x(j), z(j)));
%                         score_loo(i) = score_loo(i) + log(  KCDE_fun_loo6(x(j), z(j), j));
                        score_temp(j) = log(  KCDE_fun6(x(:,j), z(j)));
                        score_loo_temp(j) = log(  KCDE_fun_loo6(x(:,j), z(j), j));

                    end
                    score_loo(i) = mean(score_loo_temp);
                    score(i) = mean(score_temp);
                    
                    score_cv_temp = zeros(1,size(z_cv,2));
                    for j = 1:size(z_cv,2)
                        fprintf(fileID,['Approach 6 Run ',num2str(rs), ', out-of-sample: ',num2str(j/size(z_cv,2)),'\n'])
                        [~, ~,~, ~, KCDE_fun6, KCDE_fun_loo6] = find_both_h2_ultimate(x,z,z_type,'both','variable',center_x, z_cv(j),weight); 
                        score_cv_temp(j) =  log(  KCDE_fun6(x_cv(:,j), z_cv(j)));
                        

                    end
                    
                    score_cv(i) = mean(score_cv_temp);
                    
                    
                end

            end
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

    

    for i = 1:6
        if ismember(i,nums_to_do)
                for j = 1:length(zgrid)
                    
                    
                    
                    if i == 6
                        fprintf(fileID,['Approach 6 Run ',num2str(rs), ', calculate H^2: ',num2str(j/length(zgrid)),'\n'])
                        [~, ~, output6, ~, KCDE_fun6, KCDE_fun_loo6] = find_both_h2_ultimate(x,z,z_type,'both','variable',center_x, zgrid(j),weight); 
                        hz6(j) = abs(output6(2));
                        alpha6(j) = abs(output6(1));
                        f_est = @(xi) KCDE_fun6(xi, zgrid(j));
                    
                    
                    
                    
                    else
                        f_est = @(xi) KCDE_fun{i}(xi, zgrid(j));
                    end
                    
                    f_true = @(xi) fstar_true(xi,zgrid(j));
                    Hdist(i, j) = pqdistances(f_est, f_true, 2,'Hellinger', range_star(zgrid(j)));
%                     ISEdist(i, j) = pqdistances(f_est, f_true, 1,'ISE', [ min(x)-3*std(x), max(x)+3*std(x) ]);
%                     Hdist(i, j) = sqrt(Hdist(i, j) );
                end
        end
    end

    %%
%     figure;
%     for i = 1:size(hz,2)
%         if size(hz{i},1)>0
% 
%             plot( zgrid, Hdist(i,:),'linewidth',2), hold on
%         end
%     end
%     legs{1} = ['rule-of-thumb hz & hx  (',num2str(mean(Hdist(1,:)),'%.2f'),')'];
%     legs{2} = ['rule-of-thumb hz, variable hx  (',num2str(mean(Hdist(2,:)),'%.2f'),')'];
%     legs{3} = ['ML constant hz, variable hx  (',num2str(mean(Hdist(3,:)),'%.2f'),')'];
%     legs{4} = ['ML constant hz, constant hx, both optimized (',num2str(mean(Hdist(4,:)),'%.2f'),')'];
%     legs{5} = ['ML constant hz, variable hx (\alpha), both optimized  (',num2str(mean(Hdist(5,:)),'%.2f'),')'];

%     legs_used = cell(1,length(nums_to_do));
%     counter = 1;
%     for i = 1:5
%         if sum(i == nums_to_do)
%             legs_used{counter} = legs{i};
%             counter = counter + 1;
%         end
%     end
%     leg = legend(legs_used);
% 
%     set(leg,'location','southoutside','fontsize',14)
%     xlabel('$z$','interpreter','latex')
%     title('squared Hellinger distance vs $z$','interpreter','latex')
%     set(gca,'fontsize',16)


    %%
%     if strcmp(z_marginal, 'categorical')
% 
% 
% 
%         hx_vals = zeros( size(hz,2), length(zgrid));
% 
%         for i = 1:size(hz,2)
%                 if size(hz{i},1)>0
%                     for j = 1:length(zgrid)
%                         hx_vals(i,j) = hx{i}(zgrid(j));
%                     end
% 
%                 end
%         end
% 
% 
% 
%     end

%     figure;
%     for i = 1:size(hz,2)
%         if size(hz{i},1)>0
% 
%             plot( zgrid, hx_vals(i,:),'linewidth',2), hold on
%         end
%     end
%     legs{1} = ['rule-of-thumb hz & hx  (',num2str(mean(hx_vals(1,:)),'%.2f'),')'];
%     legs{2} = ['rule-of-thumb hz, variable hx  (',num2str(mean(hx_vals(2,:)),'%.2f'),')'];
%     legs{3} = ['ML constant hz, variable hx  (',num2str(mean(hx_vals(3,:)),'%.2f'),')'];
%     legs{4} = ['ML constant hz, constant hx, both optimized (',num2str(mean(hx_vals(4,:)),'%.2f'),')'];
%     legs{5} = ['ML constant hz, variable hx (\alpha), both optimized  (',num2str(mean(hx_vals(5,:)),'%.2f'),')'];
% 
%     legs_used = cell(1,length(nums_to_do));
%     counter = 1;
%     for i = 1:5
%         if sum(i == nums_to_do)
%             legs_used{counter} = legs{i};
%             counter = counter + 1;
%         end
%     end
%     leg = legend(legs_used);
% 
%     set(leg,'location','southoutside','fontsize',14)
%     xlabel('$z$','interpreter','latex')
%     title('$hx(z)$ vs $z$','interpreter','latex')
%     set(gca,'fontsize',16)

    hdists_rs(:,rs) = mean(Hdist');
    score_rs(:,rs) = score;
    score_loo_rs(:,rs) = score_loo;
    score_cv_rs(:,rs) = score_cv;
    toc


fclose(fileID);
cd results
save([fileName(1:end-4),'.mat'],'hdists_rs','score_rs','score_loo_rs','score_cv_rs','alpha6','alphas_rs')


