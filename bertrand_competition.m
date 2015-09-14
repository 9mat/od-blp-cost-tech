diary;
clear;
load goo'd solution 6.mat'

factor = 0:.02:.4;

for i = 1:length(factor);
    Data2 = Data;
    Data2.pgreal = (1+factor(i))*Data2.pgreal;
    Data2.dpm = Data2.gpm.*Data2.pgreal;
    
    Data2.Xrc(:,2) = Data2.pgreal;
    Data2.Xv(:,2) = Data2.pgreal;
    
    delta2 = delta - beta(2)*Data.pgreal + beta(2)*Data2.pgreal;
    
    Data2.c = c;
    Data2.eb = exp(eb);
    
    p_sol = zeros(size(price));
    e_sol = zeros(size(e));
    
    for t = 9:9
        fprintf(' market %d\n',t);
        index = iT == t;
        Data3.Xrc       = Data2.Xrc(index,:);
        Data3.price     = Data2.price(index,:);
        Data3.share     = Data2.share(index);
        Data3.v         = Data2.v(index,:,:);
        Data3.vprice    = Data2.vprice(index,:);
        Data3.ve        = Data2.ve(index,:);
        Data3.gpm       = Data2.gpm(index);
        Data3.dpm       = Data2.dpm(index);
        Data3.pgreal    = Data2.pgreal(index);
        Data3.iT        = Data2.iT(index);
        Data3.iF        = Data2.iF(index);
        Data3.c         = Data2.c(index);
        Data3.eb        = Data2.eb(index);
        
        [T3, ~, iT3]    = unique(Data3.iT);
        [F3, ~, iF3]    = unique([Data3.iT, Data3.iF], 'rows');
        
        Data3.iT        = iT3;
        Data3.iF        = iF3;
        
        delta3          = delta2(index);
        if i == 1
            p3              = price(index);
            e3              = e(index);
        else
            p3 = psol1;
            e3 = esol1;
        end
        
        J = length(p3);
        
        lb1 = zeros(2*J,1);
        ub1 = [1000*ones(J,1); 3*ones(J,1)];
        lb2 = zeros(J,1);
        ub2 = 1000*ones(J,1);
        
        options = optimoptions('lsqnonlin', 'Display', 'iter', 'DerivativeCheck', 'off', 'Jacobian', 'on', 'TolFun', 1e-8, 'TolX', 1e-8);
        xpe = lsqnonlin(@(x) ddpi(x(1:J), x(J+1:end), delta3, theta, Data3), [p3;e3], lb1, ub1, options);
        psol2 = lsqnonlin(@(x) ddpi_price(x, e3, delta3, theta, Data3), p3, lb2, ub2, options);
        
        %     fun  = @(x) ssq_dpi(x(1:J), x(J+1:end), delta3, theta, Data3);
        %     grad = @(x) grad_ssq_dpi(x(1:J), x(J+1:end), delta3, theta, Data3);
        %     nlcon = @(x) ddpi(x(1:J), x(J+1:end), delta3, theta, Data3);
        %     nljac = @(x) grad_ddpi(x(1:J), x(J+1:end), delta3, theta, Data3);
        %     nlrhs = zeros(size(x0));
        %     nle = zeros(size(x0));
        %
        %     optObj = optiset('display', 'iter', 'tolrfun', 1e-6, 'tolafun', 1e-6,...
        %         'maxtime', 1e5, 'solver', 'IPOPT', 'derivCheck', 'on');
        %     Opt = opti('fun', fun, 'f', grad, 'nlcon', nlcon, 'nljac', nljac, 'nlrhs', nlrhs, 'nle', nle, 'x0', x0, ...
        %         'lb', lb, 'ub', ub, 'options', optObj);
        %
        %     [x, fval] = solve(Opt);
        
        psol1 = xpe(1:J);
        esol1 = xpe(J+1:end);
    end
    
    Data31 = Data3;
    Data32 = Data3;
    
    Data31.price = psol1;
    Data31.gpm = Data31.eb - esol1;
    Data31.dpm = Data31.pgreal.*Data31.gpm;
    
    Data32.price = psol2;
    
    s1 = mean(calshare(delta3, theta, Data31), 2);
    s2 = mean(calshare(delta3, theta, Data32), 2);
    
    meangpm1(i) = sum(s1.*Data31.gpm)/sum(s1);
    meangpm2(i) = sum(s2.*Data32.gpm)/sum(s2);
    
    share1(i) = sum(s1);
    share2(i) = sum(s2);
    
    meanprice1(i) = sum(s1.*Data31.price)/sum(s1);
    meanprice2(i) = sum(s2.*Data32.price)/sum(s2);
end

%%
h=figure;
plot(factor*100, meangpm1,'b+-', factor*100, meangpm2, 'ro-')
title({'Sales-weighted GPM' 'when gasoline price increases, year 2014'}, 'FontSize', 12);
xlabel('% increase in gasoline price (relative to the actual value)');
ylabel('sales-weighted gpm');
legend( 'allowing both price change and FST adoption', 'allowing price change only');
saveTightFigure(h,'counter-factual-gpm-pgrealchange.pdf');

h=figure;
plot(factor*100, meanprice1,'b+-', factor*100, meanprice2, 'ro-')
title({'Sales-weighted prices' 'when gasoline price increases, year 2014'}, 'FontSize', 12);
xlabel('% increase in gasoline price (relative to the actual value)');
ylabel('sales-weighted price');
legend( 'allowing both price change and FST adoption', 'allowing price change only');
saveTightFigure(h,'counter-factual-price-pgrealchange.pdf');

h=figure;
plot(factor*100, share1,'b+-', factor*100, share2, 'ro-')
title({'Total market shares of new vehicles' 'when gasoline price increases, year 2014'}, 'FontSize', 12);
xlabel('% increase in gasoline price (relative to the actual value)');
ylabel('total market share');
legend( 'allowing both price change and FST adoption', 'allowing price change only','Location','northwest');
saveTightFigure(h,'counter-factual-share-pgrealchange.pdf');

diary off;