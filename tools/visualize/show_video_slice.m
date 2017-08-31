function im = show_video_slice(video, loc, gap, linecolour)
	% Function to visualize a XY, XT and YT slice of a video
	%
	% Inputs:
	% 	video: 3/4D video object
	% 	loc: Location to visualize the video slice
	% 	gap: Two tuple of pixel spacing betwene XY and XT and XY and YT images.
	% 	linecolour: Three tuple colour to draw on the location.
	% 
	% Output:
	% 	im: Image with the XY, XT and YT slices

	if ndims(video) == 3
		[H, W, T] = size(video);
		video = cat(4, video, video, video);
	else
		[H, W, T, ~] = size(video);
	end

	if isa(video, 'double')
		im = ones(H + gap(1) + T, W + gap(2) + T, 3);
	else
		im = 255*ones(H + gap(1) + T, W + gap(2) + T, 3, 'uint8');
	end

	im(1:H, 1:W, :) = squeeze(video(:, :, loc(3), :));

	linecolour = reshape(linecolour, 1, 1, 3);

	im(loc(1)-1:loc(1)+1, 1:W, :) = repmat(linecolour, 3, W, 1);
	im(1:H, loc(2)-1:loc(2)+1, :) = repmat(linecolour, H, 3, 1);
	
	im(1:H, W+gap(1)+1:end, :) = squeeze(video(:, loc(2), :, :));
	im(H+gap(2)+1:end, 1:W, :) = permute(squeeze(video(loc(1), :, :, :)), ...
										 [2 1 3]);
end
