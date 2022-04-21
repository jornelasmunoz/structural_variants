function [sparsity,threshold] = sparsity_search(vec)
    threshold = logspace(-20,1,100);
    sparsity = zeros(length(threshold),1);
    for i = 1: length(threshold)
        fhat = vec > threshold(i);
        sparsity(i) = nnz(fhat);
    end
end