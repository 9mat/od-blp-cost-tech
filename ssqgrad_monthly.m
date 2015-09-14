function [ grad ] = ssqgrad_monthly( sigma, share, X1, mpd, pg, v, iT, mark )
%SSQGRAD_MONTHLY Summary of this function goes here
%   Detailed explanation goes here
[jab, delta, ~] = jacob(sigma, share, [X1, mpd, pg], v, iT);

if any(delta > 1e30)
    grad = 1e30*ones(size(sigma));
    return
end

deltar = NaN(size(mark));
deltar(mark) = delta;

mpdr = NaN(size(mark));
mpdr(mark) = mpd;

pgr = NaN(size(mark));
pgr(mark) = pg;

jabr = NaN([numel(mark) numel(sigma)]);
jabr(mark(:), :) = jab;
jabr = reshape(jabr, [size(mark) numel(sigma)]);

Ddelta = diff(deltar, 1, 2);
Dmpd = diff(mpdr, 1, 2);
Dpg  = diff(pgr, 1, 2);

Djab = diff(jabr,1,2);
Djab = reshape(Djab, [numel(Djab)/numel(sigma) numel(sigma)]);

mark = ~isnan(Ddelta);
XX = [Dmpd(mark) Dpg(mark)];
YY = Ddelta(mark);
gamma = (XX'*XX)\(XX'*YY);

xi = YY - XX*gamma;

dd = Djab(mark,:);

grad = (2*(dd - XX/(XX'*XX)*XX'*dd)'*xi)';


end

