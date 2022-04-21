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
% =         Preparation of data variables: Simulated Data
% =========================================================================


addpath([genpath('/Users/jocelynornelasmunoz/Desktop/Research/structural_variants/'), ...
         genpath('/Users/jocelynornelas/iCloud Drive (Archive)/Desktop/UC Merced/Research/structural_variants/') ])
% -------------------------  Load Simulated Data  -------------------------


%filename = 'data/dummy_2pctNovel_2k_6n.mat';


% Varying coverage datasets for n= 10^4
filenames = ["data/10000n_5k/2Lp_4Lc/diploid_2pctNovel_60pctSim.mat",...
             "data/10000n_5k/4Lp_2Lc/diploid_2pctNovel_60pctSim.mat",...
             "data/10000n_5k/4Lp_4Lc/diploid_2pctNovel_60pctSim.mat",...
             "data/10000n_5k/4Lp_8Lc/diploid_2pctNovel_60pctSim.mat",...
             "data/10000n_5k/4Lp_16Lc/diploid_2pctNovel_60pctSim.mat",...
             "data/10000n_5k/8Lp_4Lc/diploid_2pctNovel_60pctSim.mat",...
             "data/10000n_5k/16Lp_4Lc/diploid_2pctNovel_60pctSim.mat"];
% % Varying coverage datasets for n= 10^5
%filenames = ["data/100000n_50k/2Lp_4Lc/diploid_2pctNovel_60pctSim.mat",...
% "data/100000n_50k/4Lp_2Lc/diploid_2pctNovel_60pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_2pctNovel_60pctSim.mat",...
% "data/100000n_50k/4Lp_8Lc/diploid_2pctNovel_60pctSim.mat",...
% "data/100000n_50k/4Lp_16Lc/diploid_2pctNovel_60pctSim.mat",...
% "data/100000n_50k/8Lp_4Lc/diploid_2pctNovel_60pctSim.mat",...
% "data/100000n_50k/16Lp_4Lc/diploid_2pctNovel_60pctSim.mat"]
% 
% % Varying percent novel 
% filenames = ["data/100000n_50k/4Lp_4Lc/diploid_4pctNovel_60pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_6pctNovel_60pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_8pctNovel_60pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_10pctNovel_60pctSim.mat"]
% 
% % Varying parent similarity
% filenames =                       ["data/100000n_50k/4Lp_4Lc/diploid_4pctNovel_50pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_4pctNovel_70pctSim.mat",...
% "data/100000n_50k/4Lp_4Lc/diploid_4pctNovel_80pctSim.mat]

% Haploid data
%filename = 'data/haploid_5pctNovel_10k_100n.mat';
%filename = 'data/haploid_20pctNovel_10k_100n_reproducedAPL.mat'; %reproduced data
%filename = 'lib/old/neg_binom_nov_p4_c4_5perNov.mat'; %Andrew's 5%nov 10^6n
%filename = 'lib/old/neg_binom_nov_p4_c4_20perNov.mat'; %Andrew's 20%nov
%filename = 'lib/old/dummy_data.mat';
for file = 1:length(filenames)
    filename = filenames(file);
load(filename)


if contains(filename, 'neg_binom') || contains(filename, 'dummy_data')
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
tau_vals = [0.01, 0.1, 1, 10, 100, 1000];
gamma_vals = [2, 10, 20, 100, 200, 500];
params_AUCs = zeros(length(tau_vals),length(gamma_vals), 8);
plot_flag = 1; print = 0;
if print == 1
    fprintf(['==========================================================================================\n',...
             '=                         Regularization parameters and AUCs                             =\n',...
             '= Tau        Gamma      N_total    S_total    N_parent   S_parent   N_child    S_child   =\n'])
end
for i = 1:length(tau_vals)
tau = tau_vals(i);
    for j= 1:length(gamma_vals)
        gamma = gamma_vals(j);
        reg_params =[tau; gamma];
        
        % -------------------------  Define variables  ----------------------------
        % Define true signal f, observed signal s
        s_obs = double([s_p; s_c]);
        if lower(kind) == 'haploid'
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
            
        
        elseif lower(kind) == 'diploid'
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
        miniter = 5;
        maxiter = 1000;
        tolerance = 1e-8;
        verbose = 1;
        stopcriterion = 3; 
            % 3: Relative changes in iterate
            % 4: relative changes in objective
        logepsilon = 1e-10;
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
            
            % Separate reconstruction into each signal
            if lower(kind) == 'diploid'
                fhatSPIRAL_p = 2*fhatSPIRAL(1:n)      + fhatSPIRAL(3*n+1:4*n);
                fhatSPIRAL_h = 2*fhatSPIRAL(n+1:2*n)  + fhatSPIRAL(4*n+1:5*n);
                fhatSPIRAL_n = 2*fhatSPIRAL(2*n+1:3*n)+ fhatSPIRAL(5*n+1:6*n);
                fhatSPIRAL_c = fhatSPIRAL_h + fhatSPIRAL_n;
            elseif lower(kind) == 'haploid'
                fhatSPIRAL_p = fhatSPIRAL(1:n);
                fhatSPIRAL_h = fhatSPIRAL(n+1:2*n);
                fhatSPIRAL_n = fhatSPIRAL(2*n+1:3*n);
                fhatSPIRAL_c = fhatSPIRAL_h + fhatSPIRAL_n;
            end
        
            %Separate reconstruction into each signal
            if lower(kind) == 'diploid'
                fhatNEBULA_p = 2*fhatNEBULA(1:n)      + fhatNEBULA(3*n+1:4*n);
                fhatNEBULA_h = 2*fhatNEBULA(n+1:2*n)  + fhatNEBULA(4*n+1:5*n);
                fhatNEBULA_n = 2*fhatNEBULA(2*n+1:3*n)+ fhatNEBULA(5*n+1:6*n);
                fhatNEBULA_c = fhatNEBULA_h + fhatNEBULA_n;
            elseif lower(kind) == 'haploid'
                fhatNEBULA_p = fhatNEBULA(1:n);
                fhatNEBULA_h = fhatNEBULA(n+1:2*n);
                fhatNEBULA_n = fhatNEBULA(2*n+1:3*n);
                fhatNEBULA_c = fhatNEBULA_h + fhatNEBULA_n;
            end  

            % Calculate AUCs for saving/plotting
            [X_nt,Y_nt,T_nt,AUC_nt] = perfcurve(f_true,fhatNEBULA, 1); 
            [X_st,Y_st,T_st,AUC_st] = perfcurve(f_true,fhatSPIRAL, 1);
            [X_np,Y_np,T_np,AUC_np] = perfcurve(f_p,fhatNEBULA_p, 1); 
            [X_sp,Y_sp,T_sp,AUC_sp] = perfcurve(f_p,fhatSPIRAL_p, 1);
            [X_nc,Y_nc,T_nc,AUC_nc] = perfcurve(f_c,fhatNEBULA_c, 1); 
            [X_sc,Y_sc,T_sc,AUC_sc] = perfcurve(f_c,fhatSPIRAL_c, 1);
        
            if plot_flag == 1

                % overall plot
                figure('Position', [500 500 1500 400])
                subplot(1,3,1);
                plot(X_nt,Y_nt, '-r', 'LineWidth',2); hold on
                plot(X_st, Y_st, '--b', 'LineWidth',1.5);
                legend(strcat('NEBULA = ', num2str(AUC_nt)), strcat('SPIRAL  = ', num2str(AUC_st)), 'FontSize',14, 'Location', 'southeast');
                xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
                title({'ROC Curves for total reconstruction',['\tau = ' num2str(tau),...
                    ' \gamma = ' num2str(gamma)]},'FontSize',16)
                
                
                %parent reconstruction
                subplot(1,3,2);
                plot(X_np,Y_np, '-r', 'LineWidth',2); hold on
                plot(X_sp, Y_sp, '--b', 'LineWidth',1.5);
                legend(strcat('NEBULA = ', num2str(AUC_np)), strcat('SPIRAL  = ', num2str(AUC_sp)), 'FontSize',14, 'Location', 'southeast');
                xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
                title({'ROC Curves for parent reconstruction',['\tau = ' num2str(tau),...
                    ' \gamma = ' num2str(gamma)]},'FontSize',16)
                
                %child reconstruction
                subplot(1,3,3);
                plot(X_nc,Y_nc, '-r', 'LineWidth',2); hold on
                plot(X_sc, Y_sc, '--b', 'LineWidth',1.5);
                legend(strcat('NEBULA = ', num2str(AUC_nc)), strcat('SPIRAL  = ', num2str(AUC_sc)), 'FontSize',14, 'Location', 'southeast');
                xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
                title({'ROC Curves for child reconstruction',['\tau = ' num2str(tau),...
                    ' \gamma = ' num2str(gamma)]},'FontSize',16)
            else
                %disp('The plotting flag is OFF')
            end
        params_AUCs(i, j,:) = [tau, gamma, AUC_nt, AUC_st, AUC_np, AUC_sp, AUC_nc, AUC_sc]; 
        fprintf(['-------------------------------------------------------\n',...
         '        Sparsity for vectors (nonzero counts)          \n',...
         '  NEBULA          SPIRAL          TRUTH        size\n',...
         'fhat:%6d     fhat:%6d     f:   %6d     %d\n',...
         'f_p: %6d     f_p: %6d     f_p: %6d     %d\n',...
         'f_c: %6d     f_c: %6d     f_c: %6d     %d\n',...
         'f_h: %6d     f_h: %6d     f_h: %6d     %d\n',...
         'f_n: %6d     f_n: %6d     f_n: %6d     %d\n'],...
         compute_sparsity(fhatNEBULA), compute_sparsity(fhatSPIRAL), compute_sparsity(f_true), numel(f_true),...
         compute_sparsity(fhatNEBULA_p), compute_sparsity(fhatSPIRAL_p), compute_sparsity(f_p), numel(f_p),...
         compute_sparsity(fhatNEBULA_c), compute_sparsity(fhatSPIRAL_c), compute_sparsity(f_c), numel(f_c),...
         compute_sparsity(fhatNEBULA_h), compute_sparsity(fhatSPIRAL_h), compute_sparsity(f_h), numel(f_h),...
         compute_sparsity(fhatNEBULA_n), compute_sparsity(fhatSPIRAL_n), compute_sparsity(f_n), numel(f_n))
        

        %save results
        sub_folder = char(fileparts(filename));
        save_folder = strcat('results/',sub_folder(6:end));
        if ~exist(save_folder, 'dir')
            disp('Making folder');disp(save_folder)
            mkdir(save_folder)
        end
        save_path = sprintf(strcat(save_folder,'/%stau_%sgamma_RESULTS.mat'), num2str(tau), num2str(gamma));
        save(save_path)
        if print == 1
            fprintf('= %-10.3f %-10.3f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f=\n',...
                     tau, gamma, AUC_nt, AUC_st, AUC_np, AUC_sp, AUC_nc, AUC_sc)
        end
    end
end

if print == 1
    fprintf(['=========================================================================================='])
end
format longg
params_AUCs = reshape(params_AUCs, i*j,8);
end


% Confusion Matrix
% % Not well labelled ATM, so may be hard to understand
%     figure
%     subplot(1,2,1)
%     threshold = 0.5;
%     fhatNEBULA_bin = fhatNEBULA > threshold;
%     C = confusionmat(f_true, double(fhatNEBULA_bin), 'Order', [1, 0]);
%     cm = confusionchart(C);
%     cm.Title = 'Structural Variant Confusion Matrix using NEBULA';
%     cm.RowSummary = 'row-normalized';
%     cm.ColumnSummary = 'column-normalized';
%     cm.NormalizedValues;
%     subplot(1,2,2)
%     fhatSPIRAL_bin = fhatSPIRAL > threshold;
%     C = confusionmat(f_true, double(fhatSPIRAL_bin), 'Order', [1, 0]);
%     cm = confusionchart(C);
%     cm.Title = 'Structural Variant Confusion Matrix using SPIRAL';
%     cm.RowSummary = 'row-normalized';
%     cm.ColumnSummary = 'column-normalized';
%     cm.NormalizedValues;

%end