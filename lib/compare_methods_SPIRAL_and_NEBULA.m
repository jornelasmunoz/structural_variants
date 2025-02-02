%compare_methods_SPIRAL_and_NEBULA.m
% This script compares the simulated diploid, 1 parent-1 child data using
% SPIRAL-- Poisson-based method and
% NEBULA-- Negative Binomial-based method

% This code was adapted from Melissa Spence and Mario Banuelos's SV work

clc;
clear;
close all
format longg

% =========================================================================
% =         Preprocessing of data variables: Simulated Data
% =========================================================================


addpath([genpath('/Users/jocelynornelasmunoz/Desktop/Research/structural_variants/'), ...
         genpath('/Users/jocelynornelas/iCloud Drive (Archive)/Desktop/UC Merced/Research/structural_variants/'),...
         genpath('/home/jornelasmunoz/structural_variants/')])

% -------------------------  Load Simulated Data  -------------------------

% Coverage (P,C)= (7,3), (3,7), (5,5) x erreps=0.1,0.5
% 09/26/22 Experiments
filenames = [
%            "data/100000n_5000k/3Lp_7Lc/diploid_4pctNovel_80pctSim_5e-02eps.mat"
%             "data/100000n_5000k/3Lp_7Lc/diploid_4pctNovel_80pctSim_1e-01eps.mat",
%              "data/100000n_5000k/3Lp_7Lc/diploid_4pctNovel_80pctSim_5e-01eps.mat",
%              "data/100000n_5000k/7Lp_3Lc/diploid_4pctNovel_80pctSim_1e-01eps.mat",
%              "data/100000n_5000k/7Lp_3Lc/diploid_4pctNovel_80pctSim_5e-01eps.mat",
              "data/100000n_5000k/5Lp_5Lc/diploid_4pctNovel_80pctSim_1e-01eps.mat",
%             "data/100000n_5000k/5Lp_5Lc/diploid_4pctNovel_80pctSim_5e-01eps.mat"
            ];
%filenames = ["data/old/dummy_2pctNovel_2k_6n.mat"];

% Varying coverage datasets for n= 10^4
%  filenames = "data/10000n_5k/2Lp_4Lc/diploid_2pctNovel_60pctSim.mat";%,...
%              "data/10000n_5k/4Lp_2Lc/diploid_2pctNovel_60pctSim.mat",...
%              "data/10000n_5k/4Lp_4Lc/diploid_2pctNovel_60pctSim.mat",...
%              "data/10000n_5k/4Lp_8Lc/diploid_2pctNovel_60pctSim.mat",...
%              "data/10000n_5k/4Lp_16Lc/diploid_2pctNovel_60pctSim.mat",...
%              "data/10000n_5k/8Lp_4Lc/diploid_2pctNovel_60pctSim.mat",...
%              "data/10000n_5k/16Lp_4Lc/diploid_2pctNovel_60pctSim.mat"];
% % % % Varying coverage datasets for n= 10^5
% filenames = ["data/100000n_50k/2Lp_4Lc/diploid_2pctNovel_60pctSim.mat",...
%              "data/100000n_50k/4Lp_2Lc/diploid_2pctNovel_60pctSim.mat",...
%              "data/100000n_50k/4Lp_4Lc/diploid_2pctNovel_60pctSim.mat",...
%              "data/100000n_50k/4Lp_8Lc/diploid_2pctNovel_60pctSim.mat",...
%              "data/100000n_50k/4Lp_16Lc/diploid_2pctNovel_60pctSim.mat",...
%              "data/100000n_50k/8Lp_4Lc/diploid_2pctNovel_60pctSim.mat",...
%              "data/100000n_50k/16Lp_4Lc/diploid_2pctNovel_60pctSim.mat",...
%              "data/100000n_50k/4Lp_4Lc/diploid_4pctNovel_60pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_6pctNovel_60pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_8pctNovel_60pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_10pctNovel_60pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_4pctNovel_50pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_4pctNovel_70pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_4pctNovel_80pctSim.mat"];

% % Varying percent novel 
% filenames = ["data/100000n_50k/4Lp_4Lc/diploid_4pctNovel_60pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_6pctNovel_60pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_8pctNovel_60pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_10pctNovel_60pctSim.mat"];
% 
% % Varying parent similarity
% filenames =                       ["data/100000n_50k/4Lp_4Lc/diploid_4pctNovel_50pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_4pctNovel_70pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_4pctNovel_80pctSim.mat"];

% Haploid data
%filename = 'data/haploid_5pctNovel_10k_100n.mat';
%filename = 'data/haploid_20pctNovel_10k_100n_reproducedAPL.mat'; %reproduced data
%filenames = ["lib/old/neg_binom_nov_p4_c4_5perNov.mat"]; %Andrew's 5%nov 10^6n
%filenames = ["lib/old/neg_binom_nov_p4_c4_20perNov.mat"]; %Andrew's 20%nov
%filename = 'lib/old/dummy_data.mat';
for file = 1:length(filenames)
    clearvars -except filenames file
    filename = filenames(file);
load(filename)
fprintf(strcat('\n', filename, '\n'))
                                                                                                                                                                                                                                                         
if contains(filename, 'neg_binom') || contains(filename, 'dummy_data')
    % preprocess APL data (used for code validation)
    fprintf('Using Andrews data \n')
    kind = 'haploid';
    f_p = f_p_neg_binom;
    f_h = f_c_inh_neg_binom;
    f_n = f_c_nov_neg_binom;
    f_c = f_h + f_n;
    s_p = y_p_neg_binom;
    s_c = y_c_neg_binom;
    A_c = A_c_neg_binom;
    A_p = A_p_neg_binom;
elseif contains(filename,'haploid')
    kind = 'haploid';
else
    kind= 'diploid';
end

% ---------------------  Regularization parameters  -----------------------
% Define parameters regularization parameters 
tau_vals = 1 ;%[0.01, 0.1, 1, 10, 100, 1000];
gamma_vals = 2;%[2, 10, 20, 100, 200, 500];

% initalize vector to save AUCs
% Tau        Gamma      N_total    S_total    N_parent   S_parent   N_child    S_child
params_AUCs = zeros(length(tau_vals),length(gamma_vals), 8);

% plot and print flags
plot_flag = 1; print = 0;
if print == 1
    fprintf(['==========================================================================================\n',...
             '=                         Regularization parameters and AUCs                             =\n',...
             '= Tau        Gamma      N_total    S_total    N_parent   S_parent   N_child    S_child   =\n'])
end


% -------------------- Begin Loop for Tau and Gamma -----------------------
for i = 1:length(tau_vals)
tau = tau_vals(i);
    for j= 1:length(gamma_vals)
        gamma = gamma_vals(j);
        reg_params =[tau; gamma];
        
        % -------------------------  Define variables  ----------------------------
        % Define true signal f, observed signal s
        s_obs = double([s_p; s_c]);
        if lower(kind) == "haploid"
            % Andrew order
            s_obs = double([s_c; s_p]); 
            f_true = double([f_h; f_n; f_p]);
            reg_params_all =[tau; tau*gamma; tau]; %[f_h; f_n; f_p]
            A(1:n,1:n) = A_c;
            A(1:n,n+1:2*n) = A_c;
            A(n+1:2*n, 2*n+1:3*n) = A_p;
        
            %s_obs = double([s_p; s_c]);
            %f_true = double([f_p; f_h; f_n]);
            subvectors = 3;
            N = length(f_true);
            n = N/subvectors;
            
            % Set up diagonal block matrix A (coverage matrix)
        %     A = sparse(2*n, 3*n);
        %     A(1:n,1:n) = A_p;
        %     A(n+1:(2*n),n+1:(2*n)) = A_c;
        %     A(n+1:(2*n), (2*n)+1 :(3*n)) = A_c;
            %reg_params_all corresponds to the penalty for 
            % each corresponding subvector in f_true
            %reg_params_all =[tau; tau; tau*gamma]; %[f_p; f_h; f_n]
            B = A;
            
        
        elseif lower(kind) == "diploid"
            disp('Doing diploid case')
            f_true = double([z_p; z_h; z_n; y_p; y_h; y_n]);
            subvectors = 6;
            N = length(f_true);
            n = N/subvectors; 
            
            % Set up diagonal block matrix A (coverage matrix)
            A = sparse(2*n, 6*n);
            A(1:n,1:n) = A_zp;
            A(1:n, 3*n+1:4*n) = A_yp;
            A(n+1:2*n, n+1:2*n) = A_zc;
            A(n+1:2*n, 2*n+1:3*n) = A_zc;
            A(n+1:2*n, 4*n+1:5*n) = A_yc;
            A(n+1:2*n, 5*n+1:6*n) = A_yc;
            B = A;
        
            %reg_params_all corresponds to the penalty for [z_p; z_h; z_n; y_p; y_h; y_n]
            % in that order
            reg_params_all =[tau; tau; tau*gamma; tau; tau; tau*gamma]; 
            
        %     % ALTERNATE ORDER
        %     f_true = double([z_h; z_n; z_p; y_h; y_n; y_p]);
        %     subvectors = 6;
        %     N = length(f_true);
        %     n = N/subvectors; 
        %     
        %     % Set up diagonal block matrix A (coverage matrix)
        %     A = sparse(2*n, 6*n);
        %     A(1:n,1:n)       = A_zc;
        %     A(1:n, n+1:2*n)  = A_zc;
        %     A(1:n,3*n+1:4*n) = A_yc;
        %     A(1:n,4*n+1:5*n) = A_yc;
        %     A(n+1:2*n, 2*n+1:3*n) = A_zp;
        %     A(n+1:2*n, 5*n+1:6*n) = A_yp;
        %     B = A;
        % 
        %     %reg_params_all corresponds to the penalty for [z_p; z_h; z_n; y_p; y_h; y_n]
        %     % in that order
        %     reg_params_all =[tau; tau*gamma; tau; tau; tau*gamma; tau]; 
        end
        
        % ------------------------- Prepare for algorithm -------------------------
        % Set up function handles for computing A and A^T
        AT  = @(x)A'*x;
        A   = @(x)A*x;
        
        % set maximum number of iterations, tol, and when to print to screen
        miniter = 3;
        maxiter = 70;
        tolerance = 1e-8;
        verbose = 2;
        stopcriterion = 3; 
            % 3: Relative changes in iterate
            % 4: relative changes in objective
        logepsilon = erreps; %1e-2;
        alphamin = 1e-30;
        alphamax = 1e30;
        alphaaccept = 1e30;
        
        % Simple initialization:
        % initialization of f to start with AT(s) rescaled to a 
        % least-squares fit to the mean intensity
        f_init = (sum(sum(s_obs)).*numel(AT(s_obs)))...
            ./(sum(sum(AT(s_obs))) .*sum(sum((AT(ones(size(s_obs)))))))...
            .*AT(s_obs);
        

        
        
        % =========================================================================
        % =                  SPIRAL Method reconstruction                         =
        % =========================================================================
        
            % Run SPIRAL algorithm
            [fhatSPIRAL, iterations_SPIRAL, objective_SPIRAL,...
            reconerror_SPIRAL, cputime_SPIRAL] ...
                = NEBULA(s_obs,A,reg_params, reg_params_all,'poisson',subvectors,...
                'maxiter',maxiter,...
                'Initialization',f_init,...
                'AT',AT,...
                'miniter',miniter,...
                'stopcriterion',stopcriterion,...
                'tolerance',tolerance,...
                'alphainit',1,...
                'alphamin', alphamin,...
                'alphamax', alphamax,...
                'alphaaccept',alphaaccept,...
                'logepsilon',logepsilon,...
                'saveobjective',1,...
                'savereconerror',1,...
                'savecputime',1,...
                'savesolutionpath',0,...
                'truth',f_true,...
                'verbose',verbose);
        
            
        % =========================================================================
        % =                  NEBULA Method reconstruction                         =
        % =========================================================================
        
            % Run NEBULA algorithm
            [fhatNEBULA, iterations_NEBULA, objective_NEBULA,...
            reconerror_NEBULA, cputime_NEBULA] ...
                = NEBULA(s_obs,A,reg_params, reg_params_all,'negative binomial',subvectors,...
                'maxiter',maxiter,...
                'Initialization',f_init,...
                'AT',AT,...
                'miniter',miniter,...
                'stopcriterion',stopcriterion,...
                'tolerance',tolerance,...
                'alphainit',1,...
                'alphamin', alphamin,...
                'alphamax', alphamax,...
                'alphaaccept',alphaaccept,...
                'logepsilon',logepsilon,...
                'saveobjective',1,...
                'savereconerror',1,...
                'savecputime',1,...
                'savesolutionpath',0,...
                'truth',f_true,...
                'verbose',verbose);
            
            fprintf('max fhatNEBULA = %5.5f \n max fhatSPIRAL = %5.5f \n', max(fhatNEBULA), max(fhatSPIRAL))
            % Separate SPIRAL reconstruction into each signal
            if lower(kind) == 'diploid'
                fhatSPIRAL_zp = fhatSPIRAL(1:n);
                fhatSPIRAL_zh = fhatSPIRAL(n+1:2*n);
                fhatSPIRAL_zn = fhatSPIRAL(2*n+1:3*n);
                fhatSPIRAL_yp = fhatSPIRAL(3*n+1:4*n);
                fhatSPIRAL_yh = fhatSPIRAL(4*n+1:5*n);
                fhatSPIRAL_yn = fhatSPIRAL(5*n+1:6*n);
                
            
            elseif lower(kind) == 'haploid'
                fhatSPIRAL_p = fhatSPIRAL(1:n);
                fhatSPIRAL_h = fhatSPIRAL(n+1:2*n);
                fhatSPIRAL_n = fhatSPIRAL(2*n+1:3*n);
                fhatSPIRAL_c = fhatSPIRAL_h + fhatSPIRAL_n;
            end
        
            %Separate NEBULA reconstruction into each signal
            if lower(kind) == 'diploid'
                fhatNEBULA_zp = fhatNEBULA(1:n);
                fhatNEBULA_zh = fhatNEBULA(n+1:2*n);
                fhatNEBULA_zn = fhatNEBULA(2*n+1:3*n);
                fhatNEBULA_yp = fhatNEBULA(3*n+1:4*n);
                fhatNEBULA_yh = fhatNEBULA(4*n+1:5*n);
                fhatNEBULA_yn = fhatNEBULA(5*n+1:6*n);
            
            elseif lower(kind) == 'haploid'
                fhatNEBULA_p = fhatNEBULA(1:n);
                fhatNEBULA_h = fhatNEBULA(n+1:2*n);
                fhatNEBULA_n = fhatNEBULA(2*n+1:3*n);
                fhatNEBULA_c = fhatNEBULA_h + fhatNEBULA_n;
            end  

            % Calculate AUCs for saving/plotting
            % total reconstruction
            [X_nt,Y_nt,T_nt,AUC_nt] = perfcurve(f_true,fhatNEBULA, 1); 
            [X_st,Y_st,T_st,AUC_st] = perfcurve(f_true,fhatSPIRAL, 1);
            % homozygous 
            [X_np_2,Y_np_2,T_np_2,AUC_np_2] = perfcurve(z_p,fhatNEBULA_zp, 1); 
            [X_sp_2,Y_sp_2,T_sp_2,AUC_sp_2] = perfcurve(z_p,fhatSPIRAL_zp, 1);
            [X_nc_2h,Y_nc_2h,T_nc_2h,AUC_nc_2h] = perfcurve(z_h,fhatNEBULA_zh, 1); 
            [X_sc_2h,Y_sc_2h,T_sc_2h,AUC_sc_2h] = perfcurve(z_h,fhatSPIRAL_zh, 1);
            [X_nc_2n,Y_nc_2n,T_nc_2n,AUC_nc_2n] = perfcurve(z_n,fhatNEBULA_zn, 1); 
            [X_sc_2n,Y_sc_2n,T_sc_2n,AUC_sc_2n] = perfcurve(z_n,fhatSPIRAL_zn, 1);
            
            % heterozygous 
            [X_np_1,Y_np_1,T_np_1,AUC_np_1] = perfcurve(y_p,fhatNEBULA_yp, 1); 
            [X_sp_1,Y_sp_1,T_sp_1,AUC_sp_1] = perfcurve(y_p,fhatSPIRAL_yp, 1);
            [X_nc_1h,Y_nc_1h,T_nc_1h,AUC_nc_1h] = perfcurve(y_h,fhatNEBULA_yh, 1); 
            [X_sc_1h,Y_sc_1h,T_sc_1h,AUC_sc_1h] = perfcurve(y_h,fhatSPIRAL_yh, 1);
            [X_nc_1n,Y_nc_1n,T_nc_1n,AUC_nc_1n] = perfcurve(y_n,fhatNEBULA_yn, 1); 
            [X_sc_1n,Y_sc_1n,T_sc_1n,AUC_sc_1n] = perfcurve(y_n,fhatSPIRAL_yn, 1);

            if plot_flag == 1

                % overall plot
                figure('Position', [500 500 1500 400])
                subplot(2,4,1);
                plot(X_nt,Y_nt, '-r', 'LineWidth',2); hold on
                plot(X_st, Y_st, '--b', 'LineWidth',1.5);
                legend(strcat('NEBULA = ', num2str(AUC_nt)), strcat('SPIRAL  = ', num2str(AUC_st)), 'FontSize',14, 'Location', 'southeast');
                xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
                title({'ROC Curves for total reconstruction',['\tau = ' num2str(tau),...
                    ' \gamma = ' num2str(gamma)]},'FontSize',16)
                
                
                %parent reconstruction - homozygous
                subplot(2,4,2);
                plot(X_np_2,Y_np_2, '-r', 'LineWidth',2); hold on
                plot(X_sp_2, Y_sp_2, '--b', 'LineWidth',1.5);
                legend(strcat('NEBULA = ', num2str(AUC_np_2)), strcat('SPIRAL  = ', num2str(AUC_sp_2)), 'FontSize',14, 'Location', 'southeast');
                xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
                title({'ROC Curves for parent reconstruction - homozygous',['\tau = ' num2str(tau),...
                    ' \gamma = ' num2str(gamma)]},'FontSize',16)
                
                %child reconstruction - inherited homozygous
                subplot(2,4,3);
                plot(X_nc_2h,Y_nc_2h, '-r', 'LineWidth',2); hold on
                plot(X_sc_2h, Y_sc_2h, '--b', 'LineWidth',1.5);
                legend(strcat('NEBULA = ', num2str(AUC_nc_2h)), strcat('SPIRAL  = ', num2str(AUC_sc_2h)), 'FontSize',14, 'Location', 'southeast');
                xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
                title({'ROC Curves for child reconstruction - inherited homozygous',['\tau = ' num2str(tau),...
                    ' \gamma = ' num2str(gamma)]},'FontSize',16)
                
                %child reconstruction - novel homozygous
                subplot(2,4,4);
                plot(X_nc_2n,Y_nc_2n, '-r', 'LineWidth',2); hold on
                plot(X_sc_2n, Y_sc_2n, '--b', 'LineWidth',1.5);
                legend(strcat('NEBULA = ', num2str(AUC_nc_2n)), strcat('SPIRAL  = ', num2str(AUC_sc_2n)), 'FontSize',14, 'Location', 'southeast');
                xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
                title({'ROC Curves for child reconstruction - novel homozygous',['\tau = ' num2str(tau),...
                    ' \gamma = ' num2str(gamma)]},'FontSize',16)
                
                %parent reconstruction - heterozygous
                subplot(2,4,6);
                plot(X_np_1,Y_np_1, '-r', 'LineWidth',2); hold on
                plot(X_sp_1, Y_sp_1, '--b', 'LineWidth',1.5);
                legend(strcat('NEBULA = ', num2str(AUC_np_1)), strcat('SPIRAL  = ', num2str(AUC_sp_1)), 'FontSize',14, 'Location', 'southeast');
                xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
                title({'ROC Curves for parent reconstruction - heterozygous',['\tau = ' num2str(tau),...
                    ' \gamma = ' num2str(gamma)]},'FontSize',16)
                
                %child reconstruction - inherited heterozygous
                subplot(2,4,7);
                plot(X_nc_1h,Y_nc_1h, '-r', 'LineWidth',2); hold on
                plot(X_sc_1h, Y_sc_1h, '--b', 'LineWidth',1.5);
                legend(strcat('NEBULA = ', num2str(AUC_nc_1h)), strcat('SPIRAL  = ', num2str(AUC_sc_1h)), 'FontSize',14, 'Location', 'southeast');
                xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
                title({'ROC Curves for child reconstruction - inherited heterozygous',['\tau = ' num2str(tau),...
                    ' \gamma = ' num2str(gamma)]},'FontSize',16)
                
                %child reconstruction - novel heterozygous
                subplot(2,4,8);
                plot(X_nc_1n,Y_nc_1n, '-r', 'LineWidth',2); hold on
                plot(X_sc_1n, Y_sc_1n, '--b', 'LineWidth',1.5);
                legend(strcat('NEBULA = ', num2str(AUC_nc_1n)), strcat('SPIRAL  = ', num2str(AUC_sc_1n)), 'FontSize',14, 'Location', 'southeast');
                xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
                title({'ROC Curves for child reconstruction - novel heterozygous',['\tau = ' num2str(tau),...
                    ' \gamma = ' num2str(gamma)]},'FontSize',16)
            else
                %disp('The plotting flag is OFF')
            end
%         params_AUCs(i, j,:) = [tau, gamma, AUC_nt, AUC_st, AUC_np, AUC_sp, AUC_nc, AUC_sc]; 
%         fprintf(['-------------------------------------------------------\n',...
%          '        Sparsity for vectors (nonzero counts)          \n',...
%          '  NEBULA          SPIRAL          TRUTH        size\n',...
%          'fhat:%6d     fhat:%6d     f:   %6d     %d\n',...
%          'f_p: %6d     f_p: %6d     f_p: %6d     %d\n',...
%          'f_c: %6d     f_c: %6d     f_c: %6d     %d\n',...
%          'f_h: %6d     f_h: %6d     f_h: %6d     %d\n',...
%          'f_n: %6d     f_n: %6d     f_n: %6d     %d\n'],...
%          compute_sparsity(fhatNEBULA), compute_sparsity(fhatSPIRAL), compute_sparsity(f_true), numel(f_true),...
%          compute_sparsity(fhatNEBULA_p), compute_sparsity(fhatSPIRAL_p), compute_sparsity(f_p), numel(f_p),...
%          compute_sparsity(fhatNEBULA_c), compute_sparsity(fhatSPIRAL_c), compute_sparsity(f_c), numel(f_c),...
%          compute_sparsity(fhatNEBULA_h), compute_sparsity(fhatSPIRAL_h), compute_sparsity(f_h), numel(f_h),...
%          compute_sparsity(fhatNEBULA_n), compute_sparsity(fhatSPIRAL_n), compute_sparsity(f_n), numel(f_n))
        

        %save results
        sub_folder = char(fileparts(filename));
        char_filename = char(filename);
        save_folder = strcat('results/',sub_folder(6:end));
        if ~exist(save_folder, 'dir')
            disp('Making folder');disp(save_folder)
            mkdir(save_folder)
        end
    
        save_path = sprintf(strcat(save_folder,char_filename(length(sub_folder)+1:end-4),'_%stau_%sgamma_RESULTS.mat'), num2str(tau), num2str(gamma));
        disp(save_path)
        save(save_path)
        if print == 1
            fprintf('= %-10.3f %-10.3f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f=\n',...
                     tau, gamma, AUC_nt, AUC_st, AUC_np_2, AUC_np_1, AUC_sp_2, AUC_sp_1, AUC_nc_2n,AUC_nc_2h, AUC_nc_1h,AUC_nc_1n,AUC_sc_2h,AUC_sc_2n,AUC_sc_1h,AUC_sc_1n)
        end
    end
end

if print == 1
    fprintf(['=========================================================================================='])
end
format longg
% params_AUCs = reshape(params_AUCs, length(tau_vals)*length(gamma_vals),8);
end


