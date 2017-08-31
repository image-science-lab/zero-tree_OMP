function output_video = comp_video_denoise(video, D1, k1, D0, k0, W, U, ph, pw, pt, sh, sw, st, mp_mode)
	% Function for video denoising using constrained OMP
    % approximation.
	%
	% Input:
	% 	video: Video matrix.
	% 	D1: Low level dictionary.
	% 	k1: Low level sparsity.
	% 	D0: High level dictionary.
	% 	k0: High level sparsity.
	% 	W: Downsampler.
	% 	U: Upsampler.
	% 	ph, pw, pt: Height, width, number of frames of the video patch.
	% 	sh, sw, st: Stride along the axes.
	%   mp_mode: 'lf': y = D0X0; 'hf': y = D0X0 + UD1X1
	%
	% Output:
	% 	output_video: Recovered video.
	%

	[Hv, Wv, T] = size(video);
	output_video = zeros(Hv, Wv, T);
	weights = zeros(Hv, Wv, T);

	l0_resolve = size(D0, 2)/size(D1, 2);

	for h = 1:sh:Hv - ph
		for w = 1:sw:Wv-pw
			for t = 1:st:T-pt
				idx_h = h:h+ph-1;
				idx_w = w:w+pw-1;
				idx_t = t:t+pt-1;

				% Get the video patch.
				vid_patch = video(idx_h, idx_w, idx_t);
				patch = reshape(vid_patch, ph*pw*pt, 1);
				
				% Solve for low resolution.
				[x, sup] = sparse_omp(D1, W*patch, k1, 1e-4);
				x1 = sparse(sup, ones(length(sup), 1), ones(length(sup), 1), size(D1, 2), 1);
                lf_patch = D1(:, sup)*x;
				if strcmp(mp_mode, 'lf')
					% Now solve for high resolution.
					[x, sup] = sparse_comp(D0, patch, x1, k0, l0_resolve);

					rpatch = D0(:, sup)*x;
				elseif strcmp(mp_mode, 'hf')
					[x, sup] = sparse_comp(D0, patch - U*D1*x1, x1, k0, l0_resolve);
					rpatch = D0(:, sup)*x + U*lf_patch;
				else
					fprintf('Need MP mode to be lf or hf\n');
					return;
				end
			
				weights(idx_h, idx_w, idx_t) = weights(idx_h, idx_w, idx_t) + 1;	
				output_video(idx_h, idx_w, idx_t) = ...
					output_video(idx_h, idx_w, idx_t) + reshape(rpatch, ph, pw, pt);
			end
		end
	end
	% Renormalize the video.
	output_video(weights == 0) = video(weights == 0);
	weights(weights == 0) = 1;
	output_video = output_video./weights;
end
