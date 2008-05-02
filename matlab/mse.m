function y = mse (x, ref)

y = sum (sum (sum (abs (x - ref).^2))) / sum (sum (sum (abs (ref).^2)));
