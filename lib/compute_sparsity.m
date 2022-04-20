function sparsity = compute_sparsity(f)
    %compute the sparsity of vector f
    %assumes f is a column vector
    f(f < 1e-15)=0;
    sparsity = nnz(f);
%     size_vec = size(f);
%     sparsity = size_vec(1) - nonzero_count;
end 
