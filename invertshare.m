function [delta, s] = invertshare(theta, Data)

global lastdelta

delta = lastdelta;

% contraction mapping
converged = false;
toler = 1e-10;

restart = false;
%delta = log(share) - log(outshr);

iter = 0;
while ~converged
    delta0 = delta;
    
    s = calshare(delta, theta, Data);
    delta = delta + log(Data.share) - log(mean(s,2));
    distance = max(abs(delta(:) - delta0(:)));
    
    if any((delta > 1e30) | any(isnan(delta)))
        if restart
            return
        end
        restart = true;
        delta = zeros(size(delta));
    else
        converged = distance < toler;
    end
    
    iter = iter + 1;
    if mod(iter,1000) == 0
        fprintf(' Iteration #%d, distance = %f\n', iter, distance);
    end
    
    if iter > 10000
        delta = nan(size(delta));
    end
end

lastdelta = delta;

end