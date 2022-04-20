function [tversky_index,thresholds] = tversky_index(f_pred,f_true,alpha)
thresholds = linspace(0,1);
tversky_index = zeros(length(thresholds),1);
for  i = 1:length(thresholds)
    pred = zeros(length(f_pred));
    pred = f_pred > thresholds(i);
    disp(nnz(pred))
    TP = sum(f_true.*pred);
    FN = sum(f_true.*(1-pred));
    FP = sum((1-f_true).*pred);
    tversky_index(i) = TP/ (TP + alpha*FN + (1-alpha)*FP);

end