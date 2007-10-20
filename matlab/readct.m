function v = readct (fname);

fid = fopen (fname, 'r');

sz = fread (fid, 3, 'int32').';
y = fread (fid, 2 * prod (sz), 'float32');

fclose (fid);

v = complex (y(1:2:end), y(2:2:end));
v = reshape (v, sz);
