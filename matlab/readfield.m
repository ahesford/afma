function fld = readfield (fname);

fid = fopen (fname, 'r');

len = fread (fid, 1, 'int32');
fld = fread (fid, 2 * len, 'float32');

fclose (fid);

fld = complex (fld(1:2:end), fld(2:2:end));
