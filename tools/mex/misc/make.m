sources = {'../utils.c'};
args = {'-I../', '-I.', '-largeArrayDims'};

files = {'sdot.c'};

for idx = 1:length(files)
	fprintf('Compiling %s\n', files{idx});
	mex(files{idx}, sources{:}, args{:});
end
