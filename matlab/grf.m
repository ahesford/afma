function y = grf (k, r1, r2)

R = norm (r1 - r2);

y = exp (1i * k * R) ./ R;
