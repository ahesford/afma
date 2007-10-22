function [x, y, z, c] = scatsfc (loc, field, nphi, ntheta)

[phi, theta, rho] = cart2sph (loc(:,1), loc(:,2), loc(:,3));

prange = linspace (min (phi), max (phi), nphi);
trange = linspace (min (theta), max (theta), ntheta);

[pr, tr] = meshgrid (prange, trange);

dbfield = 10 * log10 (abs (field ./ max (field)));
dbfield = dbfield - min (dbfield);

c = griddata (phi, theta, dbfield, pr, tr);
[x, y, z] = sph2cart (pr, tr, c);
