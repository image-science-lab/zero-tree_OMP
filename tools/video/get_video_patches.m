function patches = get_video_patches(video, h, w, t, npatches)
	% Function to get patch form a 3D matrix video.

	patches = zeros(h*w*t, npatches);

	[H, W, T] = size(video);
	h_indices = h:H-h; w_indices = w:W-w; t_indices = t:T-t;
	patch_idx = 1;
	for h_idx = h_indices(randperm(H-2*h))
		for w_idx = w_indices(randperm(W-2*w))
			for t_idx = t_indices(randperm(T-2*t))
				% Get a patch first.
				patch = reshape(video(h_idx-h/2:h_idx+h/2-1, ...
									  w_idx-w/2:w_idx+w/2-1, ...
									  t_idx-t/2:t_idx+t/2-1), h*w*t, 1);
				% Check its standard deviation.
				stddev = std(double(patch));
				pmean = mean(double(patch));
				if stddev >= 0.05*pmean
					% Acccept it.
					patches(:, patch_idx) = patch./norm(patch);
					patch_idx = patch_idx + 1;
				end
				if patch_idx > npatches
					return;
				end
			end
		end
	end
	patches = patches(:, 1:patch_idx-1);
end
