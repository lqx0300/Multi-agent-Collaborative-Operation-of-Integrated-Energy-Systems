function [pri_e_21,pri_e_23,Obj_MG2]=Fun_MG2(kesai_21,kesai_23,lambda_e_21,lambda_e_23,rou,price_max,price_min,P_e_21,P_e_23,C_gap)
%微网2(MG2)的分布式优化迭代模型
%% 决策变量初始化
pri_e_21=sdpvar(1,24); %微网2向微网1交互的电价
pri_e_23=sdpvar(1,24); %微网2向微网3交互的电价
%% 导入电/热负荷和电网购电电价
pri_e=[0.40*ones(1,7),0.75*ones(1,4),1.20*ones(1,3),0.75*ones(1,4),1.20*ones(1,4),0.40*ones(1,2)];

%% 约束条件
C=[];
for t=1:24
   C=[C,
      pri_e(t)>=pri_e_21(t)>=0.2,
      pri_e(t)>=pri_e_23(t)>=0.2,
     ];
end
C=[C,C_gap+sum(pri_e_21.*P_e_21+pri_e_23.*P_e_23)>=0,];
%% 目标函数
C_trade=sum(pri_e_21.*P_e_21+pri_e_23.*P_e_23);
Obj=-log(C_gap+C_trade)...
    +0.5*rou*(norm(pri_e_21-kesai_21)^2)+sum(lambda_e_21.*pri_e_21)...
    +0.5*rou*(norm(pri_e_23-kesai_23)^2)+sum(lambda_e_23.*pri_e_23);    
%% 求解器配置与求解
ops=sdpsettings('solver','mosek','verbose',0,'usex0',0);
result=solvesdp(C,Obj,ops);
%% 数据输出
if result.problem == 0     
    pri_e_21=double(pri_e_21);
    pri_e_23=double(pri_e_23);
    Obj_MG2=double(C_trade);
else
 disp('求解失败，失败原因为：')
 disp(result.info)
end


end