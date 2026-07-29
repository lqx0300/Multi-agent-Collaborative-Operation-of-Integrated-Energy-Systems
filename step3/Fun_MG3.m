function [pri_e_31,pri_e_32,Obj_MG3]=Fun_MG3(kesai_31,kesai_32,lambda_e_31,lambda_e_32,rou,price_max,price_min,P_e_31,P_e_32,C_gap)
%微网3(MG3)的分布式优化迭代模型
%% 决策变量初始化
pri_e_31=sdpvar(1,24); %微网3向微网1交互的电价
pri_e_32=sdpvar(1,24); %微网3向微网2交互的电价
%% 导入电/热负荷和电网购电电价
pri_e=[0.40*ones(1,7),0.75*ones(1,4),1.20*ones(1,3),0.75*ones(1,4),1.20*ones(1,4),0.40*ones(1,2)];

%% 约束条件
C=[];
for t=1:24
   C=[C,
      pri_e(t)>=pri_e_31(t)>=0.2,
      pri_e(t)>=pri_e_32(t)>=0.2,
     ];
end
C=[C,C_gap+sum(pri_e_31.*P_e_31+pri_e_32.*P_e_32)>=0,];
%% 目标函数
C_trade=sum(pri_e_31.*P_e_31+pri_e_32.*P_e_32);
Obj=-log(C_gap+C_trade)...
    +0.5*rou*(norm(pri_e_31-kesai_31)^2)+sum(lambda_e_31.*pri_e_31)...
    +0.5*rou*(norm(pri_e_32-kesai_32)^2)+sum(lambda_e_32.*pri_e_32);    
%% 求解器配置与求解
ops=sdpsettings('solver','mosek','verbose',0,'usex0',0);
result=solvesdp(C,Obj,ops);
%% 数据输出
if result.problem == 0     
    pri_e_31=double(pri_e_31);
    pri_e_32=double(pri_e_32);
    Obj_MG3=double(C_trade);
else
 disp('求解失败，失败原因为：')
 disp(result.info)
end


end