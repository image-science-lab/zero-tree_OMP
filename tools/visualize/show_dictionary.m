function dict_image = show_dictionary(dictionary, isrgb)
    % Function to visualize the learnt dictionary as an image.
    
    % First, convert the dictionary to patches.
    dict_patches = visualize_dict(dictionary, isrgb);
    [m, n, c, npatches] = size(dict_patches);
    
    nimg_row = ceil(sqrt(npatches));
    nimg_col = ceil(npatches/nimg_row);
    
    dict_image = zeros(nimg_row*m, nimg_col*n, c);
    
    atom_idx = 1;
    
    for r = 1:nimg_row
        for c = 1:nimg_col
            dict_image((r-1)*m+1:r*m, (c-1)*n+1:c*n, :) = ...
                dict_patches(:, :, :, atom_idx);
            atom_idx = atom_idx + 1;
			if atom_idx > npatches
				return;
			end
        end
    end
	if c == 1
		dict_image = squeeze(dict_image);
	end
end
