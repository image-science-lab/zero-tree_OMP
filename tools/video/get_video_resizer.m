function R = get_video_resizer(h_low, h_high, t_low, t_high)
	% Function to get resizer for a video vector.

	% Get an image resizer first.
	Rim = get_resizer(h_low, h_high, false);

	% Need to combine two frames. So first double up the image resizer.
    n = t_high/t_low;
	R_double = (1/n)*repmat(Rim, 1, n);

	% Now Stack it vertically.
	R = R_double;

	for idx = 1:t_low - 1
		R = blkdiag(R, R_double);
	end
end
