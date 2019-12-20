function P = pzp(type,xx,C)
% generate obserable matrices, with memoization
% type is 'z' or 'p'

load([type 'memo.mat'], xs, Ps);

for i = 1:length(xs)
	if length(xs{i}) = length(xx) && all(xs{i} == xx) && size(Ps{i}, 2) == C+1
		P = Ps{i}
		return
	end
end




end 