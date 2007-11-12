function  M = fsfrechet (grid, src, rcv, cell)

cellv = prod (cell);

rcvgrf = pi * cellv * exp (-2i * pi * rcv * grid.');

M = [];

for i = 1:size(src,1)
	srcgrf = pi * cellv * exp(-2i * pi * src(i,:) * grid.');
	grfmat = ones(size(rcvgrf,1),1) * srcgrf;
	M = [ M; grfmat .* rcvgrf ];
end
