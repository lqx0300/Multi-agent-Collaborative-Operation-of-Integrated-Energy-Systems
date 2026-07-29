function [pri_e_12,pri_e_13,Obj_MG1]=Fun_MG1(kesai_12,kesai_13,lambda_e_12,lambda_e_13,rou,price_max,price_min,P_e_12,P_e_13,C_gap)
%微网1(MG1)的分布式优化迭代模型
%% 决策变量初始化
pri_e_12=sdpvar(1,24); %微网1向微网2交互的电价
pri_e_13=sdpvar(1,24); %微网1向微网3交互的电价

%% 导入电/热负荷和电网购电电价
pri_e=[0.40*ones(1,7),0.75*ones(1,4),1.20*ones(1,3),0.75*ones(1,4),1.20*ones(1,4),0.40*ones(1,2)];

%% 约束条件
C=[];
for t=1:24
   C=[C,
      pri_e(t)>=pri_e_12(t)>=0.2
      pri_e(t)>=pri_e_13(t)>=0.2,
     ];
end
C=[C,C_gap+sum(pri_e_12.*P_e_12+pri_e_13.*P_e_13)>=0,];
%% 目标函数
C_trade=sum(pri_e_12.*P_e_12+pri_e_13.*P_e_13);
Obj=-log(C_gap+C_trade)...
    +0.5*rou*(norm(pri_e_12-kesai_12)^2)+sum(lambda_e_12.*pri_e_12)...
    +0.5*rou*(norm(pri_e_13-kesai_13)^2)+sum(lambda_e_13.*pri_e_13);    
%% 求解器配置与求解
ops=sdpsettings('solver','mosek','verbose',0,'usex0',0);
result=solvesdp(C,Obj,ops);
%% 数据输出
if result.problem == 0     
    pri_e_12=double(pri_e_12);
    pri_e_13=double(pri_e_13);
    Obj_MG1=double(C_trade);
else
 disp('求解失败，失败原因为：')
 disp(result.info)
end


end