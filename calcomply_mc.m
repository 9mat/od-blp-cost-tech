function [comply_mc, gammai] = calcomply_mc(gamma, Data)

binding = Data.comply == 0;
fined = Data.comply == -1;

gammai = zeros(size(Data.comply));
gammai(binding) = gamma;
gammai(fined) = gamma+0.055;

comply_mc = calcomply_mc_j(gammai, Data);
