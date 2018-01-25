%% sort eigen images %%

function [pc_ind, pc_sign] = sort_eigen_images(orig_V, bs_V)

	numpcs.orig = size(orig_V, 2);
	numpcs.bs   = size(bs_V  , 2);

	numpcs      = min(numpcs.orig, numpcs.bs);

	% disp(['numpcs is equal to ', num2str(numpcs)]);

	r_tmp  = corr(orig_V, bs_V) ;
	r_sign = sign(r_tmp) ;
	r_tmp  = abs( r_tmp) ; 

	% 	disp(numpcs);

	for pc = 1:numpcs

		[ii, jj] = find(r_tmp == max(r_tmp(:))) ;

		if size(ii,1) > 1
			ii = ii(1,1);
			jj = jj(1,1);
		end

		pc_ind( ii) = jj;
		pc_sign(ii) = r_sign(ii, jj) ;
		r_tmp(ii,:) = -1;
		r_tmp(:,jj) = -1;

	end

end