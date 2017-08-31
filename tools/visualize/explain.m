function imout = explain(im, params)
    % Function to write SNR, time of execution and inset a sub image on a
    % given image.
    %
    % Inputs:
    %   im: Image to write over
    %   params: Struct of parameters
    %       snrval: SNR of reconstruction. Will be written on top right
    %       timeval: Time of reconstruction. Will be written on top left
    %       topleft: Top left corner of image to inset
    %       bottomright: Bottom right corner of the image to inset
    %       zoom_factor: Amount to zoom the inset image by
    %       width: Width of the border of inset image
    %       color: Color of the border for inset image
	% 		fontsize: Size of the font.
    %
    % Output:
    %   imout: Output image with everything written.

    % I like these colours
    lightblue = [0 176 240];
    lightgreen = [146 208 80];

    % First zoom in
	if ~isempty(params.topleft)
    	imout = zoom_in(im, params.topleft, params.bottomright, ...
                    	params.zoom_factor, params.color, params.width);
	else
		imout = im;
	end

    [~, W, ~] = size(imout);

    % Now write the timing value
    if ~isempty(params.timeval)
		if params.timeval > 60
			timestring = sprintf('%.2f min', params.timeval/60);
		else
			timestring = sprintf('%.2f s', params.timeval);
		end
        imout = insertText(imout, [0 0], timestring, ...
                           'FontSize', params.fontsize, 'TextColor', 'white', ...
                           'BoxColor', lightblue, 'BoxOpacity', 1.0);
    end
    % Finally write the SNR value
    if ~isempty(params.snrval)
        imout = insertText(imout, [W 0], sprintf('%.2fdB', params.snrval), ...
                           'AnchorPoint', 'RightTop', ...
                           'FontSize', params.fontsize, 'TextColor', 'white', ...
                           'BoxColor', lightgreen, 'BoxOpacity', 1.0);
    end
end
