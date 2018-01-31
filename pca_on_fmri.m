%% TODO %%

% make sure the bootstrap ratios are all good (try to match Randy's formula)
% optimize runtime any way possible. This takes WAY too long.
% add the input SPM filenames to teh final ouput
% is ignoring the zero the right way to fix sort_eigen()? Or should I move the zero index to the end (i.e., value equal to nsub)?

% the package path

pkg_path       = '/global/home/hpc3586/JE_packages/pca_je';
func_paths     = fullfile(pkg_path);

%% user options %%

top_dir = '/global/home/hpc3586/SART_data/output/NOGO/Combined/detrend6_NOGO_sart_combined_erCVA/optimization_results/spms' ;
output  = '/global/home/hpc3586/SART_data/output_pls/detrend6_combined_clean/NOGO/pls_outcome/yng_pca.mat' ;
mask    = '/global/home/hpc3586/bin_fun_MNI152.nii.gz'

pipe = 3;
filters = {'yng', 'sNorm'} ;

% %% run the function %%

[pca_out.scans, pca_out.pca_out, pca_out.scores, pca_out.pc_var, pca_out.st_coords, pca_out.X, pca_out.Xnorm] = pca_fmri(top_dir, output, pipe, filters, mask) ;

% %% save the output %%
save(output, 'pca_fmri') ;

%%%%%%%%%%%%%%%
%% functions %%
%%%%%%%%%%%%%%%

function [runs, pca_out, scores, pc_var, st_coords, XX, XX_norm] = pca_fmri(top_dir, output, pipe, filters, mask)

	filters = ['*' strjoin(filters, '*') '*.nii'] ;

	runs = fullfile(top_dir, filters) ;
	runs = dir(runs);
	runs = {runs.name} ;

	if exist('mask')
		mask      = load_nii(mask) ;
		mask      = mask.img ;
		mask      = reshape(mask, [1, prod(size(mask))] ) ;
		st_coords = find(mask) ;

		img_info.st_coords = st_coords;
	end

	for run = 1:size(runs,2)
		run_nm  = runs{run} ;
		run_nm  = fullfile(top_dir, run_nm) ;

		run_spm = load_nii(run_nm) ;
		run_spm = run_spm.img ;
		run_spm = run_spm(:,:,:,pipe);
		img_dim = size(run_spm);

		img_info.img_dim = img_dim;

		run_spm = reshape(run_spm, [1, prod(img_dim)]) ;

		if exist('mask')
			run_spm = run_spm(st_coords);
		end

		XX(run,:) = run_spm ;

	end

	nsub = size(XX,1);

	%% main pca %%

	disp('MAIN PCA')

	% normalization

	XX_norm = zscore(XX);

	[U, Sig, V] = svd(XX_norm, 'econ') ;

	% project eigenimages onto the BOLD image subspace to obtain subject loading scores
	scores = XX_norm * V ;

	% accounted for variance
	pc_var = diag(Sig) ;
	pc_var = pc_var ./ sum(pc_var) ;

	pca_out.U   = U;
	pca_out.Sig = Sig;
	pca_out.V   = V;

	disp('DONE');

end




