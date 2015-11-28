function [margin, flag] = calmargin(s, alphai, iF)

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
    margin(index) = -Delta\share(index);
    
    [~, warningID] = lastwarn;
    if (strcmp(warningID, 'MATLAB:illConditionedMatrix') == 1) ...
            || (strcmp(warningID, 'MATLAB:singularMatrix') == 1)
        flag = false;
        return;
    end
end


end