function y = analyt (k, dx, dy, dz)

R = (3 * dx * dy * dz / (4 * pi))^(1/3);

y = (4 * pi / (k^2)) * ((1 - 1i * k * R) * exp(1i * k * R) - 1);
