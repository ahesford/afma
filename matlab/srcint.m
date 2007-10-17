function y = srcint (k, src, obs, dx, dy, dz)

a = 0.339981043584856;
b = 0.861136311594053;
wa = 0.652145154862546;
wb = 0.347854845137454;

pts = [-b -a a b];
wts = [wb wa wa wb];

y = 0;

for i = 1:4;
	for j = 1:4;
		for l = 1:4;
			offset = 0.5 * [dx*pts(i) dy*pts(j) dz*pts(l)];
			srcpt = src + offset;
			val = rcvint (k, srcpt, obs, dx, dy, dz);
			y = y + wts(i) * wts(j) * wts(l) * val;
		end
	end
end

y = y * (dx * dy * dz / 8);
