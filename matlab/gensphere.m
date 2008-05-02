function y = gensphere (filename, hi, step, obj)

maxd = hi - 0.5 * step;
d = -maxd:step:maxd;
[x, y, z] = ndgrid (d, d, d);

r = sqrt (x.^2 + y.^2 + z.^2);

y = obj * (r < hi);

writect (y, filename);
