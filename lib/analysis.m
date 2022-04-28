% Analyze results of given file
%clc; clear all; close all;
addpath([genpath("/Users/jocelynornelasmunoz/Desktop/Research/structural_variants/"), ...
         genpath("/Users/jocelynornelas/iCloud Drive (Archive)/Desktop/UC Merced/Research/structural_variants/"),...
         genpath("/home/jornelasmunoz/structural_variants/")])

% files = dir("results/10000n_5k/2Lp_4Lc/");
% idx = 7;
% result_file = files(idx+3).name;
% %result_file = "results/10000n_5k/2Lp_4Lc/0.01tau_500gamma_RESULTS.mat";
% result_file = "results/10000n_5k/2Lp_4Lc/0.1tau_500gamma_RESULTS.mat";
% result_file = "results/10000n_5k/2Lp_4Lc/1tau_500gamma_RESULTS.mat";
% result_file = "results/10000n_5k/2Lp_4Lc/10.1tau_500gamma_RESULTS.mat";
% result_file = "results/10000n_5k/2Lp_4Lc/100tau_500gamma_RESULTS.mat";
% result_file = "results/10000n_5k/2Lp_4Lc/1000tau_500gamma_RESULTS.mat";

load(result_file)
disp(result_file)

%save results
% sub_folder = char(fileparts(result_file));
% fig_folder = strcat('figs/',sub_folder(6:end));
% if ~exist(fig_folder, 'dir')
%     disp('Making folder');disp(fig_folder)
%     mkdir(fig_folder)
%end
%save_path = sprintf(strcat(fig_folder,'/%stau_%sgamma_RESULTS.mat'), num2str(tau), num2str(gamma));
%save(save_path)
%%
%-------------------------------------------------------------------------
% ---------------------  Sparsity - i.e. nonzeros ------------------------
%-------------------------------------------------------------------------

[y_St,x_St] = sparsity_search(fhatSPIRAL);
[y_Nt,x_Nt] = sparsity_search(fhatNEBULA);
[y_Sp,x_Sp] = sparsity_search(fhatSPIRAL_p);
[y_Np,x_Np] = sparsity_search(fhatNEBULA_p);
[y_Sc,x_Sc] = sparsity_search(fhatSPIRAL_c);
[y_Nc,x_Nc] = sparsity_search(fhatNEBULA_c);

[y_T,x_T] = sparsity_search(f_true);
[y_P,x_P] = sparsity_search(f_p);
[y_C,x_C] = sparsity_search(f_c);

figure('Position', [500 500 1500 400])

subplot(1,3,1);
loglog(x_T,y_T,'--'); hold on; 
loglog(x_Nt,y_Nt,'-', 'LineWidth',1); 
loglog(x_St,y_St,'-', 'LineWidth',1); hold off
%xlabel('Threshold','FontSize',16); 
ylabel('Nonzero count','FontSize',16);
title({'Total Signal Sparsity'},'FontSize',16);
legend('Truth','NEBULA', 'SPIRAL', 'FontSize',14, 'Location','southwest');
ylim([0, 10^4]);

subplot(1,3,2);
loglog(x_P,y_P,'--'); hold on; 
loglog(x_Sp,y_Sp,'-', 'LineWidth',1);
loglog(x_Np,y_Np,'-', 'LineWidth',1); hold off;
xlabel('Threshold','FontSize',16); %ylabel('Nonzero count','FontSize',16);
title({'Parent Signal Sparsity'},'FontSize',16);
legend('Truth','NEBULA', 'SPIRAL', 'FontSize',14, 'Location','southwest');
ylim([0, 10^4])

subplot(1,3,3);
loglog(x_C,y_C,'--'); hold on; 
loglog(x_Sc,y_Sc,'-', 'LineWidth',1);
loglog(x_Nc,y_Nc,'-', 'LineWidth',1); hold off;
%xlabel('Threshold','FontSize',16); ylabel('Nonzero count','FontSize',16);
title({'Child Signal Sparsity'},'FontSize',16);
legend('Truth','NEBULA', 'SPIRAL', 'FontSize',14, 'Location','southwest');
ylim([0, 10^4])

%-------------------------------------------------------------------------
%------------------   Precision-Recall Curves  ---------------------------
%-------------------------------------------------------------------------

[p_N,r_N,no_skill_N,th_N] = compute_precision_recall(f_true,fhatNEBULA);
[p_S,r_S,no_skill_S,th_S] = compute_precision_recall(f_true,fhatSPIRAL);
[p_Np,r_Np,no_skill_Np,th_Np] = compute_precision_recall(logical(f_p),fhatNEBULA_p);
[p_Sp,r_Sp,no_skill_Sp,th_Sp] = compute_precision_recall(logical(f_p),fhatSPIRAL_p);
[p_Nc,r_Nc,no_skill_Nc,th_Nc] = compute_precision_recall(logical(f_c),fhatNEBULA_c);
[p_Sc,r_Sc,no_skill_Sc,th_Sc] = compute_precision_recall(logical(f_c),fhatSPIRAL_c);

figure('Position', [500 500 1500 400])
subplot(1,3,1);
plot([r_N;0],[p_N;1],'LineWidth',1); hold on;
plot([r_S;0],[p_S;1],'LineWidth',1);
plot([0;1], [no_skill_N;no_skill_N], '--r'); hold off
title({'Total Signal Precision-Recall'},'FontSize',16);
legend('NEBULA', 'SPIRAL','No Skill', 'FontSize',14, 'Location','southwest');
ylabel('Precision','FontSize',16)

subplot(1,3,2);
plot([r_Np;0],[p_Np,;1],'LineWidth',1); hold on;
plot([r_Sp;0],[p_Sp;1],'LineWidth',1);
plot([0;1], [no_skill_Np;no_skill_Np], '--r'); hold off
title({'Parent Signal Precision-Recall'},'FontSize',16);
legend('NEBULA', 'SPIRAL','No Skill', 'FontSize',14, 'Location','southwest');
xlabel('Recall','FontSize',16)

subplot(1,3,3);
plot([r_Nc;0],[p_Nc;1],'LineWidth',1); hold on;
plot([r_Sc;0],[p_Sc;1],'LineWidth',1);
plot([0;1], [no_skill_Nc;no_skill_Nc], '--r'); hold off
title({'Child Signal Precision-Recall'},'FontSize',16);
legend('NEBULA', 'SPIRAL','No Skill', 'FontSize',14, 'Location','southwest');

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% order [tau, gamma, AUC_nt, AUC_st, AUC_np, AUC_sp, AUC_nc, AUC_sc];
params_AUCs = reshape(params_AUCs, 6*6,8);
params_AUCs = params_AUCs(params_AUCs(:,1) == tau,:);
%params_AUCs(params_AUCs ~=0);
%params_AUCs = reshape(params_AUCs, length(params_AUCs)/8,8);

figure()
scatter(params_AUCs(:,3),params_AUCs(:,4),25, 'k');hold on
scatter(params_AUCs(:,5),params_AUCs(:,6),25, 'r');
scatter(params_AUCs(:,7),params_AUCs(:,8),25, 'b');
plot([0.5,1],[0.5,1], '--'); hold off
xlabel('NEBULA'); ylabel('SPIRAL')
legend('Total signal', 'Parent signal', 'Child signal', 'Location','southeast','FontSize',14);
title('ROC AUC comparison between NEBULA and SPIRAL','FontSize',16)

figure
h = heatmap(params_AUCs(:,3:end),'FontSize',13);
h.YDisplayLabels = string(params_AUCs(:,2));
h.XDisplayLabels = ["N_{total}", "S_{total}", "N_{parent}", "S_{parent}", "N_{child}", "S_{child}"];
h.Title = strcat('ROC AUC comparisons for \tau = ', num2str(tau));