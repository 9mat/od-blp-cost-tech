function [margin, flag] = calmargin(s, alphai, iF)

warning('error', 'MATLAB:illConditionedMatrix');
warning('error', 'MATLAB:singularMatrix');
warning('error', 'MATLAB:nearlySingularMatrix');

F = max(iF);
J = length(iF);
N = size(s,2);

margin = zeros(J,1);
share = mean(s,2);

sa = bsxfun(@times, s, alphai);
meansa = mean(sa,2);

flag = true;
for f=1:F
    index = iF == f;
    Delta = diag(meansa(index)) - sa(index,:)*s(index,:)'/N;
%     singular = rcond(Delta);
%     if (singular < eps) || isnan(singular)
%         flag(f) = false;
%     end
    lastwarn('');
    
    try
        margin(index) = -Delta\share(index);
    catch
        flag = false;
        break;
    end
end

warning('on', 'MATLAB:illConditionedMatrix');
warning('on', 'MATLAB:singularMatrix');
warning('on', 'MATLAB:nearlySingularMatrix');

end