function save_dictionary_video(dict_video, filename)
    % Function to save a GIF version of the dictionary.
    f = size(dict_video, 3);
    
    for n = 1:f
      im = dict_video(:, :, n);
      [imind,cm] = gray2ind(im,256);
      if n == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end
    end
end