top_dir = '/global/home/hpc3586/SART_data/output/NOGO/Combined/detrend6_NOGO_sart_combined_erCVA/optimization_results/spms' ;
output  = '/global/home/hpc3586/SART_data/output_pls/detrend6_combined_clean/NOGO/pls_outcome/yng_testPLS.mat' ;

%% user options %%

pipe = 3;
nperm = 1000;
nboot = 1000;
filters = {'yng', 'sNorm'} ;

%% script begins %%

filters = ['*' strjoin(filters, '*') '*'] ;

runs = fullfile(top_dir, filters) ;
runs = dir(runs);
runs = {runs.name} ;

for run = 1:size(runs,2)
	run_nm  = runs{run} ;
	run_nm  = fullfile(top_dir, run_nm) ;

	run_spm = load_nii(run_nm) ;
	run_spm = run_spm.img ;
	run_spm = run_spm(:,:,:,pipe);
	img_dim = size(run_spm);
	run_spm = reshape(run_spm, [1, prod(size(run_spm))]) ;

	XX(run,:) = run_spm ;

end

%% normalization %%

% column-wise normalization

XX_norm = zscore(XX) ;

%% SVD %%

[U, Sig, V] = svd(XX_norm, 'econ') ;

% accounted for variance
pc_var = diag(Sig) ;
pc_var = pc_var ./ sum(pc_var) ;

% compute pc scores
pc_scores = U * Sig ;

%% permutation testing %%

%% bootstrap testing %%

RSalience = 0;
MSalience = 0;

for boot = 1:nboot
	
	disp([bs num2str(boot)]) ;

	isub = ceil( size(XX,1) * rand(1,size(XX,1)) );

	XX_b = XX(isub,:);

	% normalization
	XX_b = zscore(XX_b);

	% svd
	[Ub, Sigb, Vb] = svd(XX_b, 'econ') ;

	[pc_sort, pc_signs] = sort_eigen_images(V, Vb) ;

	Vb_sort = Vb(:, pc_sort);
	Vb_sort = bsxfun(@times, Vb_sort, pc_signs);

	RSalience = RSalience + Vb_sort .^ 2 ;
	MSalience = MSalience + Vb_sort      ;

end

% var(x) = E(x2) - E(x)^2
VSalience = RSalience/nboot - (MSalience/nboot) .^ 2;

% bootstrap SE
VSalience = nboot/VSalience/(nboot - 1);

% BS ratio
bs_ratio = V ./ sqrt(VSalience);

%% output %%

pls_out.dim = img_dim  ;
pls_out.X   = XX_norm  ;
pls_out.U   = U        ;
pls_out.Sig = Sig      ;
pls_out.V   = V        ;
pls_out.var = pc_var   ;
pls_out.bsr = bs_ratio ;

%% save output %%

save(output, 'pls_out');

%% sort eigen images %%

function [pc_ind, pc_sign] = sort_eigen_images(orig_V, bs_V)

	numpcs.orig = size(orig_V, 2);
	numpcs.bs   = size(bs_V  , 2);
	numpcs      = min(numpcs.orig, numpcs.bs);

	r_tmp  = corr(orig_V, bs_V) ;
	r_sign = sign(r_tmp) ;
	r_tmp  = abs( r_tmp) ; 

	for pc = 1:numpcs

		[ii, jj] = find(r_tmp == max(r_tmp(:))) ;

		pc_ind( ii) = jj;
		pc_sign(ii) = r_sign(ii, jj) ;
		r_tmp(ii,:) = -1;
		r_tmp(:,jj) = -1;

	end

end