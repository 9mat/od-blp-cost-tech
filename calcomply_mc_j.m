function [comply_mc] = calcomply_mc_j(gammai, Data)

binding = Data.comply == 0;
inactive = Data.comply == 1;
gammai(inactive) = 0;

comply_mc1 = (1-Data.gpm./Data.cagpm).*Data.cafe;
comply_mc2 = (1-Data.cagpm./Data.cagpmstd).*Data.cafestd;
comply_mc2(binding) = 0;
comply_mc = gammai.*(comply_mc1 + comply_mc2);

comply_mc(isnan(comply_mc)) = 0;
