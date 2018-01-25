%% perform PCA %%

function [Salience, pcs, ZSalience, VSalience] = run_pca(XX_norm, nboot, nsub)

	[U, Sig, V] = svd(XX_norm, 'econ') ;

	Salience = V;

	% project eigenimages onto the BOLD image subspace
	pcs = XX_norm * V ;

	% accounted for variance
	pc_var = diag(Sig) ;
	pc_var = pc_var ./ sum(pc_var) ;

	%% bootstrap testing %%

	RSalience = 0;
	MSalience = 0;

	% disp('BOOTSTRAP');
	for boot = 1:nboot
		
		isub = ceil( size(XX_norm,1) * rand(1, size(XX_norm,1)) );

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