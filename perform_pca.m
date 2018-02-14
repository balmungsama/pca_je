pipe = 2;

% ol_GO_FIX   = {'5301_run1', '8633_run2', '4996_run1', '6194_run1'};
% ol_NOGO_FIX = {'5301_run2', '6194_run1', '4450_run2', '7324_run2', '5311_run2'};

spm_ls = dir('*.nii');
spm_ls = {spm_ls(:).name}';

for ii = 1:length(spm_ls)
	spm = spm_ls{ii};	
	spm = load_nii(spm);
	spm = spm.img;
	spm = spm(:,:,:,pipe);
	spm = reshape(spm, [1, prod(size(spm))]);
	
	XX(ii,:) = spm;
end

[pca_out.coeff  , pca_out.score  , pca_out.latent  , pca_out.tsquared  , pca_out.explained  , pca_out.mu  ] = pca(XX);

pca_out.subj_ls = spm_ls;

clear spm ii

save('../NOGO_FIX_multiout_4pc.mat', 'pca_out')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outMat_path = 'C:\Users\john\OneDrive - University of Toronto\Labs\McIntosh Lab\optimization_results' ;
outMat_path = fullfile(outMat_path, 'NOGO_FIX_pca_multiout_4pc.nii') ;

template = 'C:\Users\john\Desktop\bin_fun_MNI152.nii.gz';
template = load_nii(template);

zeroMat = zeros(size(template.img));

for pc = 1:size(pca_out.coeff,2)
	zeroMat = zeroMat .* 0;
	zeroMat(:) = pca_out.coeff(:,pc);
	outMat(:,:,:,pc) = zeroMat;
end

template.img = outMat;
template.hdr.dime.dim(5) = size(pca_out.coeff,2);

save_nii(template, outMat_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groups = zeros(size(pca_out.score,1),3);
groups(1:38,3) = 1;
groups(39:end,1) = 1;

pcs = [1,2,3];



scatter(pca_out.score(:, pcs(1)), pca_out.score(:, pcs(2)), 10, groups, 'filled'); 
text( double(pca_out.score(:, pcs(1))), double(pca_out.score(:, pcs(2))), num2str([1:size(pca_out.score,1)]') )

scatter3(pca_out.score(:, pcs(1)), pca_out.score(:, pcs(2)), pca_out.score(:, pcs(3)), 10, groups, 'filled'); 

