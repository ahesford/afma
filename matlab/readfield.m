function [loc, fld] = readfield (fname);

fid = fopen (fname, 'r');

len = fread (fid, 1, 'int32');
loc = fread (fid, 3 * len, 'float32');

fld = fread (fid, 2 * len, 'float32');

fclose (fid);

loc = [loc(1:3:end) loc(2:3:end) loc(3:3:end)];

fld = complex (fld(1:2:end), fld(2:2:end));
