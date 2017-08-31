function patches = vid2mat(video, ph, pw, pt, sh, sw, st)
	% Function to convert a video to a matrix of patches.
	% Inputs:
	% 		video: 3D array of image sequence representative
	% 			of the video.
	% 		ph: Height of the patch to extract.
	% 		pw: Width of the patch to extract.
	% 		pt: Number of frames per patch.
	% 		sh: Stride along height.
	% 		sw: Stride along width.
	% 		st: Stride along time.
	%
	% Outputs:
	% 		patches: Matrix with columns representing vectorized video
	% 			patch, each of dimension ph.pw.pt x 1

	[H, W, T] = size(video);
	patches = zeros(ph*pw*pt, round(H*W*T/(sh*sw*st)));

	npatches = 1;
	for h = 1:sh:H - ph
		for w = 1:sw:W-pw
			for t = 1:st:T-pt
				idx_h = h:h+ph-1;
				idx_w = w:w+pw-1;
				idx_t = t:t+pt-1;

				vid_patch = video(idx_h, idx_w, idx_t);
				patches(:, npatches) = reshape(vid_patch, ph*pw*pt, 1);
				npatches = npatches + 1;
			end
		end
	end
	patches = patches(:, 1:npatches-1);
end
