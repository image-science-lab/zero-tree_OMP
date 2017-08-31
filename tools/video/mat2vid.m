function [video, weights] = mat2vid(patches, ph, pw, pt, sh, sw, st, H, W, T)
	% Function to convert a matrix of patches to 3D array representing
	% a video.
	%
	% Input:
	% 		patches: Matrix of patches.
	% 		ph: Height of the patch to extract.
	% 		pw: Width of the patch to extract.
	% 		pt: Number of frames per patch.
	% 		sh: Stride along height.
	% 		sw: Stride along width.
	% 		st: Stride along time.
	% 		H: Height of the video.
	% 		W: Widht of the video.
	% 		T: Number of frames of the video.
	%
	% Outputs:
	% 		video: 3D matrix representing the video.
	% 		weights: Weight matrix for normalizing the video.

	video = zeros(H, W, T);
	weights = zeros(H, W, T);

	npatches = 1;
	for h = 1:sh:H - ph
		for w = 1:sw:W-pw
			for t = 1:st:T-pt
				idx_h = h:h+ph-1;
				idx_w = w:w+pw-1;
				idx_t = t:t+pt-1;

				vid_patch = reshape(patches(:, npatches), ph, pw, pt);
				video(idx_h, idx_w, idx_t) = video(idx_h, idx_w, idx_t) ...
							+ vid_patch;
				weights(idx_h, idx_w, idx_t) = weights(idx_h, idx_w, idx_t) ...
							+ 1;
				npatches = npatches + 1;
			end
		end
	end
	weight_div = weights;
	weight_div(weight_div == 0) = 1;
	video = video./weight_div;
end
