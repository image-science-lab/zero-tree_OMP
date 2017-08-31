function coded_video = get_coded_video(video, bump_len, coding_mode)
	% Function to get coded images from a video sequence. See Hitomi et. al's
	% paper for further reference:
	% "Video from a single coded exposure photography using a learned
	%  over-complete dictionary"
	%
	% Input:
	% 	video: Matrix of video sequence.
	% 	bump_len: Gap between two samples across time axis.
	% 	coding_mode: 'random' or 'randcol'
	%
	% Output:
	% 	coded_video: Coded video. The missing values are represented using -1.

	coded_video = -ones(size(video));
	[H, W, T] = size(video);
	nbumps = floor(T/bump_len);
	chunk_len = floor(H*W/bump_len);

	for tnum = 1:nbumps
		if strcmp(coding_mode, 'random');
			all_indices = randperm(H*W);
			
			for idx = 1:bump_len
				indices = all_indices((idx-1)*chunk_len+1:idx*chunk_len);
				frame = video(:, :, (tnum-1)*bump_len + idx);
				frame_coded = -ones(H, W);
				frame_coded(indices) = frame(indices);
				coded_video(:, :, (tnum-1)*bump_len + idx) = frame_coded;
			end
		else if strcmp(coding_mode, 'randcol')
			chunk_len = floor(H/bump_len);
			all_indices = randperm(H);
			
			for idx = 1:bump_len
				indices = all_indices((idx-1)*chunk_len+1:idx*chunk_len);
				coded_video(indices, :, (tnum-1)*bump_len + idx) = ...
					video(indices, :, (tnum-1)*bump_len + idx);
			end
		end
	end
end
