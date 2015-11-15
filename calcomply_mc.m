function [comply_mc, gammai] = calcomply_mc(gamma, Data)

binding = Data.comply == 0;
fined = Data.comply == -1;

gammai = zeros(size(Data.comply));
gammai(binding) = gamma;
gammai(fined) = gamma + 0.055;

comply_mc1 = (1-Data.gpm./Data.cagpm).*Data.cafe;
comply_mc2 = (1-Data.cagpm./Data.cagpmstd).*Data.cafestd;
comply_mc2(binding) = 0;
comply_mc = gammai.*(comply_mc1 + comply_mc2);

comply_mc(isnan(comply_mc)) = 0;
