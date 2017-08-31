function create_video(video_array, video_name, frame_rate, isrgb)
	% Function to create and save a video from a sequence of images.

    if isrgb
        [~, ~, ~, t] = size(video_array);
    else
        [~, ~, t] = size(video_array);
    end
	
	output_video = VideoWriter(video_name);
	output_video.FrameRate = frame_rate;

	open(output_video);
	for idx=1:t
        if isrgb
            writeVideo(output_video, video_array(:, :, :, idx));
        else
            writeVideo(output_video, video_array(:, :, idx));
        end
	end

	close(output_video);
end
