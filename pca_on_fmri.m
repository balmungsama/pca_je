%% TODO %%

% make sure the bootstrap ratios are all good (try to match Randy's formula)
% optimize runtime any way possible. This takes WAY too long.
% add the input SPM filenames to teh final ouput
% is ignoring the zero the right way to fix sort_eigen()? Or should I move the zero index to the end (i.e., value equal to nsub)?

% the package path

pkg_path       = '/global/home/hpc3586/JE_packages/pca_je';
subj_functions = {'run_decomp.m', 'sort_eigen.m'};
func_paths     = fullfile(pkg_path, subj_functions);

%% user options %%

top_dir = '/global/home/hpc3586/SART_data/output/NOGO/Combined/detrend6_NOGO_sart_combined_erCVA/optimization_results/spms' ;
output  = '/global/home/hpc3586/SART_data/output_pls/detrend6_combined_clean/NOGO/pca_outcome/yng_parpoolPCA.mat' ;

ncpu = 4;
pipe = 3;
nboot = 100;
filters = {'yng', 'run1', 'sNorm'} ;

%% start parpool for parellel processing 

parpool(ncpu, 'AttachedFiles', func_paths);

% %% run the function %%

[pls_fmri.avg_ZSalience, pls_fmri.pca_out] = pca_fmri(top_dir, output, pipe, filters, nboot) ;

% %% save the output %%
save(output, 'pls_fmri') ;

%%%%%%%%%%%%%%%
%% functions %%
%%%%%%%%%%%%%%%

function [avg_ZSalience, pca_out, pca_main, img_dim] = pca_fmri(top_dir, output, pipe, filters, nboot)

	if ~exist('nboot')
		nboot = 1000;
	end

	filters = ['*' strjoin(filters, '*') '*.nii'] ;

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
		run_spm = reshape(run_spm, [1, prod(img_dim)]) ;

		XX(run,:) = run_spm ;

	end

	nsub = size(XX,1);

	%% leave-one-out iterations %%

	for ii = 1:nsub
		pca_loo(ii) = struct('Salience',{NaN},'pcs',{NaN}, 'ZSalience', {NaN}, ...
									'VSalience', {NaN}, 'pcs_Xo', {NaN});
	end

	parfor ii = 1:nsub

		disp(['subject ' num2str(ii)]);

		xx       = XX;
		xo       = XX(ii,:);
		xx(ii,:) = [];

		xm   = mean(xx);
		xstd = std(xx);

		% normalization
		xo = (xo - xm)./xstd;
		xx = zscore(xx);

		% running code
		[pca_loo(ii).Salience, pca_loo(ii).pcs, pca_loo(ii).ZSalience, pca_loo(ii).VSalience] = run_decomp(xx, nboot, nsub);

		pca_loo(ii).pcs_Xo = xo * pca_loo(ii).Salience ;

	end

	%% main pca %%

	disp('MAIN PCA')

	% normalization

	XX_norm = zscore(XX);

	[pca_main.Salience, pca_main.pcs, pca_main.ZSalience, pca_main.VSalience] = run_decomp(XX_norm, nboot, nsub);

	disp('DONE');
	delete(gcp('nocreate'));

	%% averaging the results across each leave-one-out iteration %%

	pca_sort = pca_loo;
	for ii = 1:nsub

		% disp(['the bs salience is = ', num2str(size(pca_loo(ii).Salience))]);

		[ind(ii, :), sg(ii, :)] = sort_eigen(pca_main.Salience, pca_loo(ii).Salience) ;

		sg_tmp  =  sg(ii, :);
		ind_tmp = ind(ii, :);
		
		sg_tmp( find(ind_tmp == 0) ) = [];
		ind_tmp(find(ind_tmp == 0) ) = [];

		pca_sort(ii).ZSalience = bsxfun(@times,pca_loo(ii).ZSalience( : , ind_tmp), sg_tmp) ;
		pca_sort(ii).Salience  = bsxfun(@times,pca_loo(ii).Salience(  : , ind_tmp), sg_tmp) ;
		pca_sort(ii).pcs       = bsxfun(@times,pca_loo(ii).pcs(       : , ind_tmp), sg_tmp) ;
		pca_sort(ii).pcs_Xo    = bsxfun(@times,pca_loo(ii).pcs_Xo(    : , ind_tmp), sg_tmp) ;

	end

	pca_out = pca_sort;

	for ii = 1:nsub
		avg_ZSalience(:,:,ii) = pca_sort.ZSalience ;
	end

	avg_ZSalience = mean(avg_ZSalience, 3) ;

end




