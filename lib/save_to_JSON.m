%save to JSON file
data.f_true = f_true;
data.fhat_SPIRAL = fhat_SPIRAL;
data.cputime_SPIRAL = cputime_SPIRAL;
data.iterations_SPIRAL = iterations_SPIRAL;
%data.objective_SPIRAL = objective_SPIRAL
data.reconerror_SPIRAL= reconerror_SPIRAL;
data.fhat_NEBULA = fhat_NEBULA;
data.cputime_NEBULA = cputime_NEBULA;
data.iterations_NEBULA = iterations_NEBULA;
%data.objective_NEBULA = objective_NEBULA
data.reconerror_NEBULA= reconerror_NEBULA;
data.f_p = f_p;
data.f_c = f_c;
data.f_h = f_h;
data.f_n = f_n;
data.s_p = s_p;
data.s_c = s_c;
data.tau = tau;
data.gamma = gamma;
data.maxiter = maxiter;
data.n = n;
data.k = k;
data.lambda_c = lambda_c;
data.lambda_p = lambda_p; 
data.pctNovel = pctNovel;
data.erreps =  erreps;
% data.pct_similarity = pct_similarity;

test_json = jsonencode(data,'PrettyPrint',true);
%test_json = jsonencode(data,'ConvertInfAndNaN',true);
filename = sprintf('reconstruction_%dsize_%dnovel.json',pctNovel*100);
fileID = fopen(filename,'w');
test_file = fprintf(fileID, test_json);

