function [ deriv, deriv2 ] = ddpi(p, e, delta, theta, Data )
%DDPI Summary of this function goes here
%   Detailed explanation goes here

Data.price = p;
Data.gpm = Data.eb - e;
Data.dpm = Data.pgreal.*Data.gpm;

s = calshare(delta, theta, Data);

alpha = theta(end-5);
lambda = theta(end-4);
sigmap = theta(end-3);
sigmae = theta(end-2);
a = theta(end-1);
b = theta(end);

margin = p - Data.c - (a*e + b*e.^2);
dcde = a + 2*b*e;

dpidp = zeros(length(p),1);
dpide = zeros(length(e),1);

if nargout > 1
    ddpidp2 = zeros(length(p));
    ddpide2 = zeros(length(p));
    ddpidpde = zeros(length(p));
    ddpidedp = zeros(length(p));
end
N = size(s,2);

dsdp = zeros(length(p));
dsde = zeros(length(e));

for t = 1:max(Data.iT)
    index = Data.iT == t;
    si = s(index,:);
    pg = mean(Data.pgreal(index));
    
    alphai = alpha*exp(sigmap*Data.vprice(index,:));
    lambdai = lambda*exp(sigmae*Data.ve(index,:));
    
    dsdp(index,index) = (diag(sum(si.*alphai,2)) - (si.*alphai)*si')/N;
    dsde(index,index) = -(diag(sum(si.*lambdai,2)) - (si.*lambdai)*si')*pg/N;
end

for f = 1:max(Data.iF)
    index = Data.iF == f;
    nf = sum(index);
    si = s(index,:);
    
    all_t = Data.iT(index);
    t = all_t(1);
    indext = Data.iT == t;
    nt = sum(indext);
    f_of_t = Data.iF(indext);
    indexf = f_of_t == f;
    
    st = s(indext,:);
    
    pg = mean(Data.pgreal(index));
    alphai = alpha*exp(sigmap*Data.vprice(index,:));
    lambdai = lambda*exp(sigmae*Data.ve(index,:));
    
    dpidp(index) = mean(si,2) + dsdp(index,index)*margin(index);
    dpide(index) = -dcde(index).*mean(si,2) + dsde(index,index)*margin(index);
    
    if nargout > 1
        sij = permute(si, [1 3 4 2]);
        sik = permute(si, [3 1 4 2]);
        sir = permute(st, [3 4 1 2]);
        
        eyejk = eye(nf);
        eyejr = permute(eye(nf), [1 3 2]);
        eyekr = permute(eye(nf), [3 1 2]);
        eyejkr = zeros(nf,nf,nf);
        eyejkr(1:(1+nf+nf*nf):end) = 1;
        
        Ijr = zeros(nf,1,nt);
        Ikr = zeros(1,nf,nt);
        Ijkr = zeros(nf,nf,nt);
        
        Ijr(:,:,indexf) = eyejr;
        Ikr(:,:,indexf) = eyekr;
        Ijkr(:,:,indexf) = eyejkr;
        
        temp1 = bsxfun(@times, sir, eyejk);
        temp2 = bsxfun(@times, sik, Ijr);
        temp3 = bsxfun(@times, sik, Ikr);
        
        temp4 = 2*bsxfun(@times, sik, sir);
        temp5 = bsxfun(@plus, Ijkr, temp4  - temp3);
        temp5 = bsxfun(@minus, temp5, temp1 + temp2);
        ddsidmu2 = bsxfun(@times, sij, temp5);
        
        ddsdp2 = mean(bsxfun(@times, ddsidmu2, permute(alphai.^2, [1 3 4 2])),4);
        ddsde2 = mean(bsxfun(@times, ddsidmu2, permute(lambdai.^2, [1 3 4 2]))*pg^2,4);
        ddsdpde = -mean(bsxfun(@times, ddsidmu2, permute(alphai.*lambdai, [1 3 4 2]))*pg,4);
        
        ddpidp2(index,indext) = ...
            shiftdim(sum(bsxfun(@times, margin(index), ddsdp2),1),1) + dsdp(index, indext);
        ddpidp2(index,index) = ddpidp2(index,index) + dsdp(index,index)';
        
        ddpide2(index,indext) = ...
            shiftdim(sum(bsxfun(@times, margin(index), ddsde2),1),1) ...
            - bsxfun(@times, dcde(index), dsde(index,indext));
        ddpide2(index,index) = ddpide2(index,index) ...
            - diag(mean(si,2))*b*2 ...
            - bsxfun(@times, dcde(index)', dsde(index,index)');
        
        temp = shiftdim(sum(bsxfun(@times, margin(index), ddsdpde),1),1);
        ddpidpde(index,indext) = temp + dsde(index,indext);
        ddpidpde(index,index) = ddpidpde(index,index) - bsxfun(@times, dcde(index)', dsdp(index,index)');
        
        ddpidedp(index,indext) = temp - bsxfun(@times, dcde(index), dsdp(index,indext));
        ddpidedp(index,index) = ddpidedp(index,index) + dsde(index,index)';
    end
end

deriv = [dpidp; dpide]*1e4;

if nargout > 1
    deriv2 = [ddpidp2 ddpidpde; ddpidedp ddpide2]*1e4;
end

end

