function im_dict = visualize_dict(dictionary, isrgb)
    % Function to visualize the learnt dictionary.
    
    [m, n] = size(dictionary);
    
	if isrgb
		c = 3
	else
		c = 1
	end
	patch_dim = sqrt(m/c);
    im_dict = zeros(patch_dim, patch_dim, c, n);
    
    for idx = 1:1:n
        patch = reshape(dictionary(:,idx), patch_dim, patch_dim, c);
        patch = 255.0*(patch - min(patch(:)))/(max(patch(:)) - min(patch(:)));
        im_dict(:,:,:,idx) = patch;
    end
end
