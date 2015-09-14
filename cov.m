function [ V ] = cov( theta, beta, Data )
%COV Summary of this function goes here
%   Detailed explanation goes here

    function u = residuals(theta, beta)
        alpha = theta(end-5);
        lambda = theta(end-4);
        sigmap = theta(end-3);
        sigmae = theta(end-2);
        a = theta(end-1);
        b = theta(end);
        
        % invert for delta
        [delta, s] = invertshare(theta, Data);
        
        iF = Data.iF;
        c = zeros(size(iF));
        N = size(s, 2);
        alphai = alpha*exp(sigmap*Data.vprice);
        
        for f=1:max(iF)
            index = iF == f;
            si = s(index, :);
            ss = Data.share(index);
            sv = si.*alphai(index,:);
            Delta  = (diag(sum(sv,2)) - sv*si')/N;
            c(index) = Data.price(index) + Delta\ss;
        end
        
        
        part1 = zeros(size(c));
        part2 = zeros(size(c));
        
        lambdai = lambda*exp(sigmae*Data.ve);
        margin = Data.price - c;
        for f=1:max(iF)
            index = iF == f;
            si = s(index, :);
            part1(index) = Data.pgreal(index).*margin(index,:).*sum(si.*lambdai(index,:),2)/N;
            part2(index) = (lambdai(index,:).*si)*(si'*margin(index)).*Data.pgreal(index)/N;
        end
        
        c_e = -(part1-part2)./Data.share;
        e = (c_e-a)/(2*b);
        eb = log(Data.gpm + e);
        c = c - (a*e + b*e.^2);
        y = [delta; log(c); eb];
        
        u = y - Data.X*beta;
    end

% finite difference
for k = 1:length(theta)
    theta1 = theta;
    theta2 = theta;
    
    epsilon = theta(k)*1e-6;
    theta1(k) = theta1(k) - epsilon;
    theta2(k) = theta2(k) + epsilon;
    
    u1 = residuals(theta1, beta);
    u2 = residuals(theta2, beta);
    
    du(:,k) = (u2-u1)/(2*epsilon);
end

% du/dbeta = X
% du(:,end-1)=[];
D = real([du Data.X]);

u = residuals(theta, beta);
J = length(Data.price);
S = zeros(size(Data.Z,2));
for j = 1:J
    uj = u(j:J:end);
    Zj = Data.Z(j:J:end,:);
    ujZj = sum(bsxfun(@times,uj,Zj),1);
    S = S + ujZj'*ujZj;
end
S = S/J;


DZ = D'*Data.Z;
DZZZ = DZ/Data.ZZ;
DZZZZD = DZZZ*DZ';
A =DZZZZD\DZZZ;
V = J*A*S*A';
end

