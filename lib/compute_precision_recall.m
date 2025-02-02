function [precision,recall, no_skill,thresholds] = compute_precision_recall(f_true,f_pred)
% ------------------------------------------------------------------------------
% Description adapted from Python sklearn.metrics.precision_recall_curve
%     Compute precision-recall pairs for different probability thresholds.
%     Note: this implementation is restricted to the binary classification task.
%
%     Let TP = true positives, FP = false positives, FN = false negatives, 
%     then we define 
%       precision = TP/(TP+FP) and recall = TP/(TP+FN)
%     Intuitively, precision is the ability of the classifier not to label as 
%     positive a sample that is negative. The recall is the ability of the 
%     classifier to find all the positive samples.
%
%
%   INPUTS:
%       f_true : nx1 vector where n = number of samples
%       f_pred : nx1 vector of probability estimates of the positive class
%       thresholds : vector of thresholds provided by perfcurve or manually
%                    defined
%
%   OUTPUTS:
%       precision :  (n_thresholds+1)x1 vector s.t. element i is the
%                       precision of predictions with score >= thresholds[i]
%       recall :     (n_thresholds+1)x1 vector s.t. element i is the
%                       recall of predictions with score >= thresholds[i]
%       thresholds : (n_thresholds)x1 vector. Increasing thresholds on the decision 
%                       function used to compute precision and recall. 
%                       n_thresholds <= length(unique(f_pred)).

if (length(f_true) ~= length(f_pred))
    error('scores and target must have same length');
end
    total_pos = sum(f_true);
    total_neg = length(f_true)-sum(f_true);
    no_skill = length(f_true(f_true==1))/length(f_true);

    thresholds = unique(f_pred);
    precision = zeros(length(thresholds),1);
    recall = zeros(length(thresholds),1);

    for i = 1: length(thresholds)
        fhat = f_pred >= thresholds(i);
        TP = sum(f_true(fhat));
        FP = (1-f_true);
        FP = sum(FP(fhat)); 
        precision(i) = TP/(TP+FP);
        recall(i) = TP/total_pos;

    end
end