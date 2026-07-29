clc
clear
U_no=[31746,74783,42627];
C_co=[16115,65280,40843];
gap = U_no-C_co;
trade = sdpvar(1,3);
f = -(gap(1)+trade(1))*(gap(2)+trade(2))*(gap(3)+trade(3));
c = [];
c = [c,
    gap(1)+trade(1)>=0;
    gap(2)+trade(2)>=0;
    gap(3)+trade(3)>=0;
    trade(1)+trade(2)+trade(3)==0]
ops = sdpsettings('solver','mosek','verbose',0)

result = optimize(c, f, ops)

if result.problem == 0     
    trade = value(trade)
    v = gap+trade
else
 disp('求解失败，失败原因为：')
 disp(result.info)
end
