%% TODO %%

% make sure the bootstrap ratios are all good (try to match Randy's formula)
% optimize runtime any way possible. This takes WAY too long.
% add the input SPM filenames to teh final ouput

%% user options %%

top_dir = '/global/home/hpc3586/SART_data/output/NOGO/Combined/detrend6_NOGO_sart_combined_erCVA/optimization_results/spms' ;
output  = '/global/home/hpc3586/SART_data/output_pls/detrend6_combined_clean/NOGO/pls_outcome/yng_testPLS.mat' ;

pipe = 3;
nboot = 5;
filters = {'yng', 'sNorm'} ;

%% run the function %%

[pls_fmri.avg_ZSalience, pls_fmri.pls_out] = pca_fmri(top_dir, output, pipe, filters, nboot) ;

%% save the output %%
save(output, 'pls_fmri') ;

%%%%%%%%%%%%%%%
%% functions %%
%%%%%%%%%%%%%%%

function [avg_ZSalience, pls_out] = pca_fmri(top_dir, output, pipe, filters, nboot)

	if ~exist('nboot')
		nboot = 1000;
	end

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

	%% leave-one-out iterations %%

	for ii = 1:size(XX,1)
		disp(['SUBJECT ', num2str(ii)]);

		xx       = XX;
		xo       = xx(ii,:);
		xx(ii,:) = [];

		xm   = mean(xx);
		xstd = std(xx);

		% normalization
		xo = (xo - xm)./xstd;
		xx = zscore(xx);

		% running code
		[pls_loo(ii).Salience, pls_loo(ii).pcs, pls_loo(ii).ZSalience, pls_loo(ii).VSalience] = run_pca(xx, nboot);

		pls_loo(ii).pcs_Xo = xo * pls_loo(ii).Salience ;

	end

	%% main pca %%

	disp('MAIN PCA')

	% normalization

	XX_norm = zscore(XX);

	[pls_main.Salience, pls_main.pcs, pls_main.ZSalience, pls_main.VSalience] = run_pca(XX_norm, nboot);

	%% averaging the results across each leave-one-out iteration %%

	pls_sort = pls_loo;
	for ii = 1:size(XX,1)

		% disp(['pls_main: ' num2str(size(pls_main.Salience   ))]);
		% disp([' pls_loo: ' num2str(size(pls_loo(ii).Salience))]);
		disp(['ii = ' num2str(ii)]);

		[ind(ii, :), sg(ii, :)] = sort_eigen_images(pls_main.Salience, pls_loo(ii).Salience) ;
		
		disp(['size(pls_loo)               = ', num2str(size(pls_loo))]);
		disp(['size(pls_loo(ii).ZSalience) = ', num2str(size(pls_loo(ii).ZSalience))]);
		disp(['size(pls_loo(ii).Salience ) = ', num2str(size(pls_loo(ii).Salience ))]);
		disp(['size(pls_loo(ii).pcs      ) = ', num2str(size(pls_loo(ii).pcs      ))]);
		disp(['size(pls_loo(ii).pcs_Xo   ) = ', num2str(size(pls_loo(ii).pcs_Xo   ))]);
		disp(['                  ind(ii,:) = ', num2str(ind(ii, : ))]);

		signs = sg(ii ,:);
		index = ind(ii,:);
		signs(find(index == 0)) = [];
		index(find(index == 0)) = [];

		disp(['                  index     = ', num2str(index)]);
		disp(['                  signs     = ', num2str(signs)]);

		% pls_sort(ii).ZSalience = bsxfun(@times,pls_loo(ii).ZSalience( : , ind(ii,: )), sg(ii,: )) ;
		% pls_sort(ii).Salience  = bsxfun(@times,pls_loo(ii).Salience(  : , ind(ii,: )), sg(ii,: )) ;
		% pls_sort(ii).pcs       = bsxfun(@times,pls_loo(ii).pcs(       : , ind(ii,: )), sg(ii,: )) ;
		% pls_sort(ii).pcs_Xo    = bsxfun(@times,pls_loo(ii).pcs_Xo(    : , ind(ii,: )), sg(ii,: )) ;

		pls_sort(ii).ZSalience = bsxfun(@times,pls_loo(ii).ZSalience( : , index), signs) ;
		pls_sort(ii).Salience  = bsxfun(@times,pls_loo(ii).Salience(  : , index), signs) ;
		pls_sort(ii).pcs       = bsxfun(@times,pls_loo(ii).pcs(       : , index), signs) ;
		pls_sort(ii).pcs_Xo    = bsxfun(@times,pls_loo(ii).pcs_Xo(    : , index), signs) ;

	end

	pls_out = pls_sort;

	for ii = 1:size(XX,1)
		avg_ZSalience_X(:,:,ii) = pls_sort.ZSalience ;
	end

	avg_ZSalience = mean(avg_ZSalience, 3) ;

end

%% perform PCA %%

function [Salience, pcs, ZSalience, VSalience] = run_pca(XX_norm, nboot)

	[U, Sig, V] = svd(XX_norm, 'econ') ;

	Salience = V;

	% project eigenimages onto the BOLD image subspace
	pcs = XX_norm * V ;

	% accounted for variance
	pc_var = diag(Sig) ;
	pc_var = pc_var ./ sum(pc_var) ;

	% % compute pc scores
	% pc_scores = U * Sig ;

	%% permutation testing %%

	%% bootstrap testing %%

	RSalience = 0;
	MSalience = 0;

	% disp('BOOTSTRAP');
	for boot = 1:nboot
		
		disp(['	bs ' num2str(boot)]) ;

		isub = ceil( size(XX_norm,1) * rand(1,size(XX_norm,1)) );

		XX_b = XX_norm(isub,:);

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
	VSalience = nboot * VSalience/(nboot - 1);

	% BS ratio
	ZSalience = V ./ sqrt(VSalience);
end

%% sort eigen images %%

function [pc_ind, pc_sign] = sort_eigen_images(orig_V, bs_V)

	numpcs.orig = size(orig_V, 2);
	numpcs.bs   = size(bs_V  , 2);
	numpcs      = min(numpcs.orig, numpcs.bs);

	disp(['numpcs is equal to ', num2str(numpcs)]);

	r_tmp  = corr(orig_V, bs_V) ;
	r_sign = sign(r_tmp) ;
	r_tmp  = abs( r_tmp) ; 

	for pc = 1:numpcs

		[ii, jj] = find(r_tmp == max(r_tmp(:))) ;

		% disp(['[ii, jj] = ', num2str([ii, jj])]);

		if ii > size(r_sign,1)
			disp('ii is too large');
			disp(['    ii = ' num2str(size(r_sign,  1))]);
			disp(['r_sign = ', num2str(size(r_sign, 1))]);
		end

		if jj > size(r_sign,2)
			disp('ii is too large');
			disp(['    ii = ' num2str(size(r_sign,  2))]);
			disp(['r_sign = ', num2str(size(r_sign, 2))]);
		end

		% disp(['the relevant r_sign is ', num2str(r_sign(ii, jj))]);

		pc_ind( ii) = jj;
		pc_sign(ii) = r_sign(ii, jj) ;
		r_tmp(ii,:) = -1;
		r_tmp(:,jj) = -1;

	end

end