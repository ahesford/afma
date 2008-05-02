function sliceplot (v, slices, clim)

px = size(v(:,:,1));
axlim = [1 px(1) 1 px(2)];

nslice = length(slices);
ncol = floor (sqrt (nslice));
nrow = ceil (nslice / ncol);

for i = 1:nslice
	subplot (sprintf('%d%d%d', nrow, ncol, i));
	slice (v, [], [], slices(i));
	title (sprintf('z = %d pixels', slices(i)));
	axis (axlim);
	shading flat;
	set (gca, 'CLim', clim);
	colorbar vert;
end

for i = 1:nrow
	j = (i - 1) * ncol + 1;
	subplot (sprintf('%d%d%d', nrow, ncol, j));
	ylabel ('y (pixels)');
end

for i = 1:ncol
	j = (nrow - 1) * ncol + i;
	subplot (sprintf('%d%d%d', nrow, ncol, j));
	xlabel ('x (pixels)');
end
