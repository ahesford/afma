function y = writect (v, fname);

sz = size (v);
ln = prod (sz);

v = reshape (v, ln, 1);
y = zeros (2 * ln, 1);

y(1:2:end) = real(v);
y(2:2:end) = imag(v);

fid = fopen (fname, 'w');

fwrite (fid, sz, 'int32');
fwrite (fid, y, 'float32');

fclose (fid);
