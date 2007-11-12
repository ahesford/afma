function y = buildgrid (nx, ny, nz, bmin, cell)

n = nx * ny * nz;

y = [];

for gi = 0:n-1
	i = mod (gi, nx);
	j = mod (floor(gi / nx), ny);
	k = floor (gi /(nx * ny));

	cen = bmin + cell .* ([i j k] + 0.5);
	y = [ y; cen ];
end
