function [coded_video, code_block] = code_video(video, block_dim)
    % Function to get random projection of a video to create coded images.
    % Instead of completely random coding, a block is coded instead, and 
    % this block is repeated over the complete video volume.
    % Input:
    %   video: 3D matrix of video frames
    %   block_dim: [h, w, t]: dimensions of the block to code.
    %  
    % Output:
    %   coded_video: 3D matrix with coding. Missing pixels are represented
    %       by -1
    %   code_block: 3D matrix block used for coding the complete video.
    
    [H, W, T] = size(video);
    h = block_dim(1); w = block_dim(2); t = block_dim(3);
    
    coded_video = -ones(H, W, T);
    code_block = zeros(h, w, t);
    
    % Now create the code block, such that, along the third dimension,
    % there is only one 1.
    
    % We don't want the most efficient implentation. Let's just go with a
    % naive one.
    for h_idx = 1:h
        for w_idx = 1:w
            % Pick an index from 1 to t
            t_idx = randi([1, t]);
            code_block(h_idx, w_idx, t_idx) = 1;
        end
    end
    
    % Now repeat the block to cover the video.
    coding_block = repmat(code_block, ceil(H/h), ceil(W/w), ceil(T/t));
    coding_block = coding_block(1:H, 1:W, 1:T);
    
    % Now code the video.
    coded_video(coding_block == 1) = video(coding_block == 1);
end