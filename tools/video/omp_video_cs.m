function output_video = omp_video_cs(video, D, k, ph, pw, pt, sh, sw, st)
	% Function for video compressive sensing using OMP approximation.
	%
	% Input:
	% 	video: Video matrix. The missing entries should be marked by -1
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
				
				% Get the sampler matrix.
				missing_idx = find(patch < 0);
				available_idx = find(patch >= 0);
				A = sparse(1:length(available_idx), available_idx, ...
			   	ones(length(available_idx),1), length(available_idx), length(patch));

				% Use the generic patch solver now.
				recov_patch = omp_patch_solver(A, patch, D, k);
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
