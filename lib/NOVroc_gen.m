function [ FPRv, TPRv ] = NOVroc_gen( fhat_recon, f_true, thresh )

% NOVroc_gen takes in SPIRAL reconstruction, true signal, and thresholding
% values and outputs the false and true positive rates to build ROC curves.

% Initialize vectors to store True positives and False positives values
T_pv = zeros(length(thresh),1);
F_pv = zeros(length(thresh),1);

T_nv = zeros(length(thresh),1);
F_nv = zeros(length(thresh),1);

TPRv = zeros(length(thresh),1);
FPRv = zeros(length(thresh),1);

% Determine size of true and reconstructed signals
n = length(fhat_recon);

for i = 1:length(thresh)
    
    f_thresh = fhat_recon >= thresh(i);
   
    
    T_p = 0;
    F_p = 0;
    
    T_n = 0;
    F_n = 0;   
    
    for j = 1:n
        
    if f_true(j)==1 && f_thresh(j)==1 %&& fhatSPIRAL_uc(j)~=0
        T_p = T_p + 1;
        
    elseif f_true(j)==0 && f_thresh(j)==1 %&&fhatSPIRAL_uc(j)~=0
        F_p = F_p + 1;
        
    elseif  f_true(j)==0 && f_thresh(j)==0 %&& fhatSPIRAL_uc(j)==0
        T_n = T_n + 1;   
        
    elseif f_true(j)==1 && f_thresh(j)==0 %&& fhatSPIRAL_uc(j)==0 
        F_n = F_n + 1;    
    end
    
    end
    
    T_pv(i) = T_p;
    F_pv(i) = F_p;
    
    T_nv(i) = T_n;
    F_nv(i) = F_n;
    
    
    TPRv(i) = T_p/(T_p + F_n);
    FPRv(i) = F_p/(F_p + T_n);
    
end

end

