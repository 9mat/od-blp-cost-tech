function [ ssq ] = ssq_monthly(sigma, share, X1, mpd, pg, v, iT, mark)
%SSQ_MONTHLY Summary of this function goes here
%   Detailed explanation goes here

delta = invertshare(sigma, share, [X1 mpd pg], v, iT);

deltar = NaN(size(mark));
deltar(mark) = delta;

mpdr = NaN(size(mark));
mpdr(mark) = mpd;

pgr = NaN(size(mark));
pgr(mark) = pg;

Ddelta = diff(deltar, 1, 2);
Dmpd = diff(mpdr, 1, 2);
Dpg  = diff(pgr, 1, 2);

mark = ~isnan(Ddelta);
XX = [Dmpd(mark) Dpg(mark)];
YY = Ddelta(mark);
gamma = (XX'*XX)\(XX'*YY);

ssq = sum((YY-XX*gamma).^2);
end

