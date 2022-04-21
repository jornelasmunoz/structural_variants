function [precision,recall, thresholds] = compute_precision_recall(f_true,f_pred)
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

    thresholds = unique(f_pred);
    precision = zeros(length(threshold),1);
    for i = 1: length(thresholds)
        fhat = vec > thresholds(i);
        
    end
end