function [comply_mc, gammai] = calcomply_mc(gamma, Data)

binding = Data.comply == 0;
fined = Data.comply == -1;

gammai = zeros(size(Data.comply));
gammai(binding) = gamma(1);
gammai(fined) = gamma(end);

comply_mc = calcomply_mc_j(gammai, Data);
