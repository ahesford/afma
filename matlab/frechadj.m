function y = frechadj (mag, fld, grid, locs, cell)

rds = grid * locs.';

y = prod(cell) * pi * exp(-2i * pi * rds) * mag;
y = y .* fld;
