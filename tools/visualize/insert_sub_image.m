function fig_handle = insert_sub_image(im, imsub, x1, y1, x2, y2, color)
	% Insert the sub image in the image in the bottom right corner.

    if length(size(im)) == 3
        [xsub, ysub, ~] = size(imsub);
        [x, y, ~] = size(im);

        % First insert the image.
        imout = im;
        imout(x-xsub+1:x, y-ysub+1:y, :) = imsub;
    else
        [xsub, ysub] = size(imsub);
        [x, y] = size(im);

        % First insert the image.
        imout = im;
        imout(x-xsub+1:x, y-ysub+1:y) = imsub;
        
        % Replicate to form an RGB image
        imout = cat(3, imout, imout, imout);
    end
    
	% Insert a rectangle around the sub image and the sub image region
	% in the real image.
	fig_handle = figure('visible', 'off');
	imshow(imout);

	line([y-ysub, y], [x-xsub, x-xsub], 'color', 'r', 'LineWidth', 3);
	line([y-ysub, y-ysub], [x-xsub, x], 'color', 'r', 'LineWidth', 3);
	line([y-ysub, y], [x, x], 'color', 'r', 'LineWidth', 3);
	line([y, y], [x-xsub, x], 'color', 'r', 'LineWidth', 3);

	line([y1, y1], [x1, x2], 'color', 'r', 'LineWidth', 2);
	line([y1, y2], [x1, x1], 'color', 'r', 'LineWidth', 2);
	line([y1, y2], [x2, x2], 'color', 'r', 'LineWidth', 2);
	line([y2, y2], [x1, x2], 'color', 'r', 'LineWidth', 2);

end
