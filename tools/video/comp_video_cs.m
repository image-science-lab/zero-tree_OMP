function output_video = comp_video_cs(video, D1, k1, D0, k0, W, U, ph, pw, pt, sh, sw, st, mp_mode)
	% Function for video compressive sensing using constrained OMP
    % approximation.
	%
	% Input:
	% 	video: Video matrix. The missing entries should be marked by -1
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

	[H, W, T] = size(video);
	output_video = zeros(H, W, T);
	weights = zeros(H, W, T);

	for h = 1:sh:H - ph
		for w = 1:sw:W-pw
			for t = 1:st:T-pt
				idx_h = h:h+ph-1;
				idx_w = w:w+pw-1;
				idx_t = t:t+pt-1;

				% Get the video patch.
				vid_patch = video(idx_h, idx_w, idx_t);
				patch = reshape(vid_patch, ph*pw*pt, 1);
				
				% Get the sampler matrix.
				missing_idx = find(patch < 0);
                available_idx = find(patch >= 0);
                row_idx = 1:length(available_idx);
                A = sparse(row_idx, available_idx, ones(length(available_idx),1), length(available_idx), length(patch));

				% Use the generic patch solver now.
                if strcmp(mp_mode, 'lf')
                    recov_patch = comp_patch_solver(A, patch, D1, k1, D0, k0, W, U, mp_mode);
                else
                    recov_patch = rcomp(A, patch, D1, k1, D0, k0, W, U, 1);
                end
				weights(idx_h, idx_w, idx_t) = weights(idx_h, idx_w, idx_t) + 1;
				% Wherever we know the values, keep it intact.
				%patch(missing_idx) = recov_patch(missing_idx);
				output_video(idx_h, idx_w, idx_t) = ...
					output_video(idx_h, idx_w, idx_t) + reshape(recov_patch, ph, pw, pt);
			end
		end
	end
	% Renormalize the video.
	output_video(weights == 0) = video(weights == 0);
	weights(weights == 0) = 1;
	output_video = output_video./weights;
end
