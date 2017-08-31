function output_video = omp_video_denoise(video, D, k, ph, pw, pt, sh, sw, st)
	% Function for video denoising using OMP approximation.
	%
	% Input:
	% 	video: Video matrix.
	% 	D: Dictionary.
	% 	k: Sparsity.
	% 	ph, pw, pt: Height, width, number of frames of the video patch.
	% 	sh, sw, st: Stride along the axes.
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

				% Use OMP to recover patch.
				[x, sup] = sparse_omp(D, patch, k, 1e-4);
				weights(idx_h, idx_w, idx_t) = weights(idx_h, idx_w, idx_t) + 1;
				% Wherever we know the values, keep it intact.
				patch = D(:, sup)*x;
				output_video(idx_h, idx_w, idx_t) = ...
					output_video(idx_h, idx_w, idx_t) + reshape(patch, ph, pw, pt);
			end
		end
	end
	% Renormalize the video.
	output_video(weights == 0) = video(weights == 0);
	weights(weights == 0) = 1;
	output_video = output_video./weights;
end
