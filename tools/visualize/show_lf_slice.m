function im = show_lf_slice(lf, loc, linecolour, scale)
	% Function to visualize the central sub aperture view and epipolar slices
	% for a given location.
	%
	% Inputs:
	% 	lf: 4/5D light field.
	% 	loc: Location to visualize epipolar slices from.
	% 	linecolour: Three tuple colour to draw on the location.
	% 	scale: Scale to expand epipolar slice by
	% 
	% Output:
	% 	im: Image with central subaperture view, epipolar slices and lines drawn
	% 		to show location of epipolar slices.

	if ndims(lf) == 4
		[U, V, H, W] = size(lf);
		lf = cat(5, lf, lf, lf);
	else
		[U, V, H, W, ~] = size(lf);
	end

	if isa(lf, 'double')
		im = ones(H + (1+scale)*U, W + (1+scale)*V, 3);
	else
		im = 255*ones(H + (1+scale)*U, W + (1+scale)*V, 3, 'uint8');
	end

	im(1:H, 1:W, :) = squeeze(lf(round((U+1)/2), round((V+1)/2), :, :, :));

	linecolour = reshape(linecolour, 1, 1, 3);

	im(loc(1)-1:loc(1)+1, 1:W, :) = repmat(linecolour, 3, W, 1);
	im(1:H, loc(2)-1:loc(2)+1, :) = repmat(linecolour, H, 3, 1);
	
	epi_im = permute(squeeze(lf(:, 1, :, loc(2), :)), [2 1 3]);
	im(1:H, W+V+1:end, :) = imresize(epi_im, [H, scale*U]);

	epi_im = squeeze(lf(1, :, loc(1), :, :));
	im(H+U+1:end, 1:W, :) = imresize(epi_im, [scale*V, W]);
end
