function [ V ] = cov( theta, beta, Data )
%COV Summary of this function goes here
%   Detailed explanation goes here

    function u = residuals(theta, beta)
        params = getParams(theta);
        [delta, s] = invertshare(theta, Data);
        
        % marginal compliance cost
        binding = Data.comply == 0;
        fined = Data.comply == -1;
        
        gammai = zeros(size(Data.comply));
        gammai(binding) = params.gamma;
        gammai(fined) = params.gamma + 0.05;
        
        comply_mc1 = (1-Data.gpm./Data.cagpm).*Data.cafe;
        comply_mc2 = (1-Data.cagpm./Data.cagpmstd).*Data.cafestd;
        comply_mc2(binding) = 0;
        comply_mc = gammai.*(comply_mc1 + comply_mc2);
        
        comply_mc(isnan(comply_mc)) = 0;
                
        alphai = params.alpha*exp(params.sigmap*Data.vprice);
        margin = calmargin(s, alphai, Data.iF);
        c = Data.price - margin - comply_mc;
        logc = log(c);
        logc(c<=0) = 1e-30;
        y = [delta; logc];
        u = y - Data.X*beta;
    end

% finite difference
du = zeros(length(Data.price)*2, numel(theta));
for k = 1:length(theta)
    theta1 = theta;
    theta2 = theta;
    
    epsilon = theta(k)*1e-7;
    theta1(k) = theta1(k) - epsilon;
    theta2(k) = theta2(k) + epsilon;
    
    u1 = residuals(theta1, beta);
    u2 = residuals(theta2, beta);
    
    du(:,k) = (u2-u1)/(2*epsilon);
end

% du/dbeta = X
% du(:,end-1)=[];
D = [du Data.X];

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

