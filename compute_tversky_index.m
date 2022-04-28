function [tversky, thresholds] = compute_tversky_index(f_true,f_pred, alpha)

if (length(f_true) ~= length(f_pred))
    error('scores and target must have same length');
end
    total_pos = sum(f_true);
    total_neg = length(f_true)-sum(f_true);

    thresholds = unique(f_pred);
    tversky = zeros(length(thresholds),1);

    for i = 1: length(thresholds)
        fhat = f_pred >= thresholds(i);
        TP = sum(f_true(fhat));
        FP = sum(1-f_true(fhat)); 
        FN = sum(1-f_pred(logical(f_true)));
        tversky(i) = TP/(TP + alpha*FN + (1-alpha)*FP);
        
    end
end