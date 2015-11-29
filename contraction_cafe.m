function [ p, margin, s, gammaj, cafe, i, distance ] = ...
    contraction_cafe( theta, delta, c, Data, gammaj, p0, settings )
%CONTRACTION_CAFE Summary of this function goes here
%   Detailed explanation goes here

maxiter = settings.maxitercafe;
toler = settings.tolercafe;

% stepsize = 0.1;
% iter = 100;

% default input
% if nargin < 8; maxiter = 500; end;
% if nargin < 7; stepsize = 0.1; end;
if nargin < 6; p0 = Data.price; end;

fleet = Data.fleet;
cafestd = Data.cafestd;
cafestdf = cafestd(Data.fleetind);
p = p0;

minmpg = accumarray(fleet, 1./Data.gpm*100, [], @min);
maxmpg = accumarray(fleet, 1./Data.gpm*100, [], @max);
complyf = Data.comply(Data.fleetind);
notbinding = (minmpg > cafestdf) | (maxmpg < cafestdf);
maybinding = ~(notbinding(fleet) | (Data.comply == -1));

complyf( (complyf==0) & (minmpg > cafestdf)) = 1;
complyf( (complyf==0) & (maxmpg < cafestdf)) = -1;
Data.comply = complyf(fleet);
binding = Data.comply == 0;
inactive = Data.comply == 1;

gammaj(minmpg(fleet) > cafestdf(fleet)) = 0;
gammaj(maxmpg(fleet) < cafestdf(fleet)) = mean(gammaj(Data.comply==-1));

gammaj0 = gammaj;

    function [p, direction, distance, cafef, maxstep, iter, distance_bertrand] = f(gammaj, p)
        index = any(isnan(p));
        p(:,index) = p0(:,index);
        [p, margin, s, iter, ~, distance_bertrand ] ...
            = contraction_bertrand(theta, delta, c, Data, gammaj, p, settings);
        share = mean(s,2);
        cafef = accumarray(fleet, share)./accumarray(fleet, share.*(Data.gpm/100));
        
        direction = cafef - cafestdf;
        index = binding & (gammaj > 0) & (direction(fleet) > 0);
        maxstep = max(direction(fleet(index))./gammaj(index));
        maxstep = min(maxstep, 0.1);
        
        binding2 = maybinding; % | (inactive & (direction(fleet) < 0));
        
        activebinding = binding2 & ((gammaj > 0) | (direction(fleet) <= 0));
        distance = max(abs([0;direction(fleet(activebinding))]));
    end

% [p, direction, distance, ~, maxstep, iter, distance_bertrand] = f(gammaj, p, 10000);
% fprintf('**** Starting: distance = %f, step = %f, bertrand iter = % 5d, distance = %f\n', distance, maxstep, iter, distance_bertrand);
% 

stepsize = @(r,v) -norm(r)/norm(v);
% fprintf('     - starting distance = %f, bertrand iter = % 5d, distance = %f \n', distance, iter, distance_bertrand);
tic;
for i = 1:maxiter    
    step = 0.1;
%     fprintf('**** CAFE iter #%5d\n', i);

    [p, direction, distance, cafef, maxstep, iter, distance_bertrand] = f(gammaj, p);
    gammaj2 = gammaj;
    binding2 = maybinding; % | (inactive & (direction(fleet) < 0));
    gammaj2(binding2) = gammaj2(binding2) - step*direction(fleet(binding2));
    gammaj2(gammaj2 < 0) = 0;

    r = gammaj2 - gammaj;
    
    distance = max(abs(r)/step);
    fprintf('CAFE iter #%4d, dist = %f, bertrand iter = % 4d, dist = %f, time = %.1fs \n', i, distance, iter, distance_bertrand, toc);
    if distance < toler; break; end;
    
    [p2, direction, distance3, ~, maxstep, iter3, distance_bertrand3] = f(gammaj2, p);
    gammaj3 = gammaj2;
    binding2 = maybinding; % | (inactive & (direction(fleet) < 0));
    gammaj3(binding2) = gammaj3(binding2) - step*direction(fleet(binding2));
    gammaj3(gammaj3 < 0) = 0;
    
    if (distance_bertrand3 < 1e-3) && (distance_bertrand < 1e-3)
        v = (gammaj3 - gammaj2) - r;
        a = stepsize(r,v);
        gammaj = gammaj - 2*a*r + a^2*v;
%         distance = min(distance, distance3);
        p = p2;
        distance = max(abs(gammaj3 - gammaj2)/step);
%         fprintf('     - 2nd step, distance = %f, bertrand iter = % 4d, dist = %f, time = %.1f secs \n', distance, iter3, distance_bertrand3);
       
    else
        gammaj = (gammaj + gammaj0)/2;
%         if distance_bertrand < 1e-3
%             gammaj = gammaj + r*0.1;
%         else
%             gammaj(maybinding) = gammaj(maybinding).*(1 + (rand(sum(maybinding),1)-0.5)*0.1);
%         end
    end
    gammaj(gammaj<0) = 0;
%     [p, direction, distance, cafef, maxstep, iter, distance_bertrand] = f(gammaj, p);

%     if i == 6; keyboard; end;
%     fprintf('     - 2nd step, distance = %f, bertrand iter = % 5d, distance = %f \n', distance, iter, distance_bertrand);
   
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

% fprintf('**     Done CAFE iteration % 4d, distance = %f\n', i, distance);
% fprintf('************************************************************************\n');
end

