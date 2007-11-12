function y = frechet (dc, fld, grid, locs, cell)

crt = dc .* fld;

rds = locs * grid.';

y = prod(cell) * pi * exp(2i * pi * rds) * crt;
