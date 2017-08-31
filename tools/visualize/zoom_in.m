function zoomed_im = zoom_in(im, topleft, bottomright, zoom_factor, color, width)
	% Function to zoom in on a tiny patch in the image and paste it in
	% the bottom right corner of the image.
	%
	% Inputs:
	% 	im: Image. Can be grayscale or RGB.
	% 	topleft: Top-left coordinates of the patch to extract.
	% 	bottomright: Bottom-right coordinates of the patch to extract.
	% 	zoom_factor: The amount by which to zoom the image.
	% 	color: Color of the border around patch and the zoomed in image.
    %   width: Width of the line to draw.
	%
	% Outputs:
	% 	zoomed_im: Output image.

	% Extract the patch first.
	x1 = topleft(1); y1 = topleft(2);
	x2 = bottomright(1); y2 = bottomright(2);

	if ndims(im) == 2
		patch = im(x1:x2, y1:y2);
	else
		patch = im(x1:x2, y1:y2, :);
	end
	imsub = imresize(patch, zoom_factor);

	% Convert image to an RGB copy.
	if ndims(im) == 2
		imsub = imsub(:, :, [1 1 1]);
		im = im(:, :, [1 1 1]);
	end
	
	% Get the shape of the sub image.
	[xsub, ysub, ~] = size(imsub);

	zoomed_im = im;
	zoomed_im(end-xsub+1:end, end-ysub+1:end, :) = imsub;
    k = width;
	% Now draw lines around patch.
	zoomed_im(x1:x2, y1-k:y1+k, 1) = color(1);
	zoomed_im(x1:x2, y1-k:y1+k, 2) = color(2);
	zoomed_im(x1:x2, y1-k:y1+k, 3) = color(3);

	zoomed_im(x1:x2, y2-k:y2+k, 1) = color(1);
	zoomed_im(x1:x2, y2-k:y2+k, 2) = color(2);
	zoomed_im(x1:x2, y2-k:y2+k, 3) = color(3);

	zoomed_im(x1-k:x1+k, y1:y2, 1) = color(1);
	zoomed_im(x1-k:x1+k, y1:y2, 2) = color(2);
	zoomed_im(x1-k:x1+k, y1:y2, 3) = color(3);
	
	zoomed_im(x2-k:x2+k, y1:y2, 1) = color(1);
	zoomed_im(x2-k:x2+k, y1:y2, 2) = color(2);
	zoomed_im(x2-k:x2+k, y1:y2, 3) = color(3);

	% Draw around zoomed patch.
	zoomed_im(end-xsub+1:end, end-ysub-k:end-ysub+k, 1) = color(1);
	zoomed_im(end-xsub+1:end, end-ysub-k:end-ysub+k, 2) = color(2);
	zoomed_im(end-xsub+1:end, end-ysub-k:end-ysub+k, 3) = color(3);

	zoomed_im(end-xsub-k:end-xsub+k, end-ysub+1:end, 1) = color(1);
	zoomed_im(end-xsub-k:end-xsub+k, end-ysub+1:end, 2) = color(2);
	zoomed_im(end-xsub-k:end-xsub+k, end-ysub+1:end, 3) = color(3);

	zoomed_im(end-xsub+1:end, end-2:end, 1) = color(1);
	zoomed_im(end-xsub+1:end, end-2:end, 2) = color(2);
	zoomed_im(end-xsub+1:end, end-2:end, 3) = color(3);

	zoomed_im(end-2:end, end-ysub+1:end, 1) = color(1);
	zoomed_im(end-2:end, end-ysub+1:end, 2) = color(2);
	zoomed_im(end-2:end, end-ysub+1:end, 3) = color(3);
end
