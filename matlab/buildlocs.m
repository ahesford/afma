function locs = buildlocs (rad, tmin, tmax, ntheta, pmin, pmax, nphi)

theta = linspace (tmin * pi / 180, tmax * pi / 180, ntheta);
phi = linspace (pmin * pi / 180, pmax * pi / 180, nphi);

x = rad * cos(phi.') * sin(theta);
y = rad * sin(phi.') * sin(theta);
z = rad * ones(size(phi.')) * cos(theta);

x = reshape (x, ntheta * nphi, 1);
y = reshape (y, ntheta * nphi, 1);
z = reshape (z, ntheta * nphi, 1);

locs = [x y z];
