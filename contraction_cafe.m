function [ p, margin, s, gammaj, cafe, i, distance ] = ...
    contraction_cafe( theta, delta, c, Data, gammaj, p0, stepsize, max_iter )
%CONTRACTION_CAFE Summary of this function goes here
%   Detailed explanation goes here

% stepsize = 0.1;
% iter = 100;

% default input
if nargin < 8; max_iter = 500; end;
if nargin < 7; stepsize = 0.1; end;
if nargin < 6; p0 = Data.price; end;

fleet = Data.fleet;
cafestd = Data.cafestd;
cafestdf = cafestd(Data.fleetind);
p = p0;

minmpg = accumarray(fleet, 1./Data.gpm*100, [], @min);
maxmpg = accumarray(fleet, 1./Data.gpm*100, [], @max);
complyf = Data.comply(Data.fleetind);
complyf( (complyf==0) & (minmpg > cafestdf)) = 1;
complyf( (complyf==0) & (maxmpg < cafestdf)) = -1;
Data.comply = complyf(fleet);
binding = Data.comply == 0;

gammaj(Data.comply == 1) = 0;

toler = 1e-5;

    function [p, direction, distance, cafef, maxstep, iter, distance_bertrand] = f(gammaj, p, maxiter)
        index = any(isnan(p));
        p(:,index) = p0(:,index);
        [p, margin, s, iter, ~, distance_bertrand ] ...
            = contraction_bertrand(theta, delta, c, Data, gammaj, p, maxiter);
        share = mean(s,2);
        cafef = accumarray(fleet, share)./accumarray(fleet, share.*(Data.gpm/100));
        
        direction = cafef - cafestdf;
        index = binding & (gammaj > 0) & (direction(fleet) > 0);
        maxstep = max(direction(fleet(index))./gammaj(index));
        maxstep = min(maxstep, 0.1);
        
        activebinding = binding & ((gammaj > 0) | (direction(fleet) <= 0));
        distance = max(abs([0;direction(fleet(activebinding))]));
    end

% [p, direction, distance, ~, maxstep, iter, distance_bertrand] = f(gammaj, p, 10000);
% fprintf('**** Starting: distance = %f, step = %f, bertrand iter = % 5d, distance = %f\n', distance, maxstep, iter, distance_bertrand);
% 

stepsize = @(r,v) -norm(r)/norm(v);
[p, direction, distance, cafef, maxstep, iter, distance_bertrand] = f(gammaj, p, 10000);
fprintf('     - starting distance = %f, bertrand iter = % 5d, distance = %f \n', distance, iter, distance_bertrand);

for i = 1:max_iter    
    step = 0.1;
    fprintf('**** CAFE iter #%5d\n', i);
    

    gammaj2 = gammaj;
    gammaj2(binding) = gammaj2(binding) - step*direction(fleet(binding));
    gammaj2(gammaj2 < 0) = 0;

    r = gammaj2 - gammaj;
    
    distance = max(abs(r)/step);
    if distance < toler; break; end;
    
    [p2, direction, distance3, ~, maxstep, iter3, distance_bertrand3] = f(gammaj2, p, 2000);
    fprintf('     - 1st step, distance = %f, bertrand iter = % 5d, distance = %f \n', distance3, iter3, distance_bertrand3);
    gammaj3 = gammaj2;
    gammaj3(binding) = gammaj3(binding) - step*direction(fleet(binding));
    gammaj3(gammaj3 < 0) = 0;
    
    if distance_bertrand < 1e-3
        v = (gammaj3 - gammaj2) - r;
        a = stepsize(r,v);
        gammaj = gammaj - 2*a*r + a^2*v;
        distance = min(distance, distance3);
        p = p2;
    else
        gammaj = gammaj + r;
    end
    gammaj(gammaj<0) = 0;
    [p, direction, distance, cafef, maxstep, iter, distance_bertrand] = f(gammaj, p, 10000);

%     if i == 6; keyboard; end;
    fprintf('     - 2nd step, distance = %f, bertrand iter = % 5d, distance = %f \n', distance, iter, distance_bertrand);
   
%     if (distance2 > lastdistance); break; end
    
%     lastdistance = distance2;
%     lastgammaj = gammaj2;
%     step = step*1.05;

%     while true
%         gammaj2 = gammaj;
%         gammaj2(binding) = gammaj2(binding) - step*direction(fleet(binding));
%         gammaj2(gammaj2 < 0) = 0;
%         [p, direction, distance2, cafef, maxstep, iter, distance_bertrand] = f(gammaj2, p, 10000);
%         
%         fprintf('   reduced step = %f, distance = %f, bertrand iter = % 5d, distance = %f \n', step, distance2, iter, distance_bertrand);
%         if (distance2 > lastdistance); break; end
%         
%         lastdistance = distance2;
%         lastgammaj = gammaj2;
%         step = step*1.05;
%     end
    
%     gammaj = lastgammaj;
%     distance = lastdistance;
    
%     gammaj(gammaj<0) = 0;
%     if distance < toler; break; end
%     
%     if olddistance - distance < 0
%         stepsize = stepsize*0.9;
%         gammaj = oldgammaj;
%         gammaj(binding) = gammaj(binding) - stepsize*error(fleet(binding));
%         gammaj(gammaj<0) = 0;
%         %         if mod(i, 100) == 0;
%         fprintf('****** ###### CAFE iteration % 4d, reduced stepsize = %f, old distance = %f, distance = %f\n', i, stepsize, olddistance, distance);
%         %         end
%     else
%         oldgammaj = gammaj;
%         olddistance = distance;
%         p = p2;
%         gammaj(binding) = gammaj(binding) - stepsize*error(fleet(binding));
%         gammaj(gammaj<0) = 0;
%         stepsize = stepsize*1.05;
%         %         if mod(i, 100) == 0;
%         fprintf('****** CAFE iteration % 4d, bertrand #iter = % 4d, error = %f, stepsize = %f, old distance = %f, distance = %f\n', i, iter, stepsize, distance_bertrand, olddistance, distance);
%         %         end
%     end
end

cafe = cafef(fleet);

fprintf('**     Done CAFE iteration % 4d, distance = %f\n', i, distance);
fprintf('************************************************************************\n');
end
