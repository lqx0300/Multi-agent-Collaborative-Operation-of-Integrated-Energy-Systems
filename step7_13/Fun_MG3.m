function [P_e_31,P_e_32,Obj_MG3,p3,h3,c3]=Fun_MG3(z_31,z_32,lambda_e_31,lambda_e_32,rou)
%微网3(MG3)的分布式优化迭代模型
%% 决策变量初始化
L_e=sdpvar(1,24); %微网经过需求响应后实际的电负荷
L_h=sdpvar(1,24); %微网经过需求响应后实际的热负荷
P_e_cut=sdpvar(1,24);   %微网的可削减电负荷
P_e_tran=sdpvar(1,24);  %微网的可转移电负荷
P_h_DR=sdpvar(1,24);    %微网的可削减热负荷
E_bat=sdpvar(1,24);     %微网中的储电设备的储电余量
P_batc=sdpvar(1,24);    %储电设备的充电功率
P_batd=sdpvar(1,24);    %储电设备的放电功率
U_abs=binvar(1,24);     %储电设备的放电状态位,取1时为放电,0为未放电
U_relea=binvar(1,24);   %储电设备的充电状态位,取1时为充电,0为未充电
P_e_pv=sdpvar(1,24); %风力的实际出力值
P_e_GT=sdpvar(1,24); %燃气轮机的发电功率
P_h_GT=sdpvar(1,24); %燃气轮机的产热功率
P_h_GB=sdpvar(1,24); %余热锅炉的产热功率
P_buy=sdpvar(1,24);  %微网向外电网的购买的电功率
P_sell=sdpvar(1,24); %微网向外电网的售出的电功率
Gas_GT=sdpvar(1,24); %GT的耗气量
Gas_GB=sdpvar(1,24); %GB的耗气量
Gas=sdpvar(1,24);    %系统的总耗气量
P_e_31=sdpvar(1,24);
P_e_32=sdpvar(1,24);
p3=zeros(9,24);
h3=zeros(3,24);
c3=zeros(1,6);
%% 导入电/热负荷和电网购电电价
L_e0=[1180,1073,1196,1165,1165,1165,1196,1503,1733,2147,2193,2407,2469,2653,2852,2208,2285,2300,3220,3220,2423,1993,1779,1518]*1.5; 
L_h0=[1494,1448,1325,1302,1317,1494,1594,1833,2080,2211,2311,2388,2303,2526,2434,2326,2164,2126,2118,2519,1841,1625,1494,1394]*1.5;
Predict_pv=[0,0,0,0,0,0,663,1084,1903,2277,2386,2480,2402,2168,2012,1474,998,0,0,0,0,0,0,0]*2;
pri_e=[0.40*ones(1,7),0.75*ones(1,4),1.20*ones(1,3),0.75*ones(1,4),1.20*ones(1,4),0.40*ones(1,2)];
grid_sw=[0.2*ones(1,24)]; 
%% 约束条件
C=[];
%微网的电/热负荷需求响应部分
for t=1:24
    C=[C,
       L_e(t)==L_e0(t)+P_e_cut(t)+P_e_tran(t), %微网的电负荷功率平衡约束
       L_h(t)==L_h0(t)-P_h_DR(t), %微网的热负荷功率平衡约束
       -0.15*L_e0(t)<=P_e_cut(t)<=0, %微网的可削减电功率上下限约束
       -0.15*L_e0(t)<=P_e_tran(t)<=0.15*L_e0(t), %微网的可转移电功率上下限约束
       0<=P_h_DR(t)<=0.2*L_h0(t), %微网的可削减热功率上下限约束
      ];
end
C=[C,sum(P_e_tran)==0,]; %转移的电负荷总量为0约束
%微网的储电设备约束部分
%储能电站荷电状态连续性约束
C=[C,E_bat(1)==800+0.95*P_batc(1)-P_batd(1)/0.96,]; %1时段约束
for t=2:24
    C=[C,E_bat(t)==E_bat(t-1)+0.95*P_batc(t)-P_batd(t)/0.96,]; %储电设备容量变化约束
end
%储能容量大小约束
for t=1:24
    C=[C,500<=E_bat(t)<=1800,];  %储电量上下限约束
end
%始末状态守恒
C=[C,E_bat(24)==800,];
%储能电站的充放电功率约束,Big-M法进行线性化处理
M=800; %这里的M是个很大的数
for t=1:24
    C=[C,
       0<=P_batc(t)<=500,
       0<=P_batc(t)<=U_abs(t)*M,     
       0<=P_batd(t)<=600,      
       0<=P_batd(t)<=U_relea(t)*M,
       U_abs(t)+U_relea(t)<=1,
      ];
end
%CHP机组以及GB等设备运行约束
for t=1:24
    C=[C,
       P_h_GT(t)==(1-0.35)/0.35*0.83*P_e_GT(t), %燃气轮机热电联产功率约束
       0<=P_e_GT(t)<=3000, %燃气轮机发电功率上下限约束
       0<=P_h_GB(t)<=500, %余热锅炉产热功率上下限约束
       0<=P_e_pv(t)<=Predict_pv(t), %风力发电上下限约束
      ];
end
%微网的热/电负荷平衡约束
for t=1:24
    C=[C,
       P_e_GT(t)+P_e_pv(t)+P_buy(t)+P_batd(t)==P_batc(t)+L_e(t)+P_e_31(t)+P_e_32(t)+P_sell(t),
       P_h_GT(t)+P_h_GB(t)==L_h(t),
      ];
end
%变量非负性等约束
for t=1:24
    C=[C,
       P_buy(t)>=0,
       0<=P_sell(t)<=2000,
       2000>=P_e_31(t)>=-2000,
       0>=P_e_32(t)>=0,
      ];
end
%
for t=1:24
    C=[C,
       P_e_GT(t)==0.35*9.7*Gas_GT(t), %GT耗气量约束
       P_h_GB(t)==0.9*9.7*Gas_GB(t), %GB耗气量约束
       Gas(t)==Gas_GT(t)+Gas_GB(t),  %总耗气量约束
      ];
end
E_co2=0.55*sum(P_e_GT)+0.65*sum(P_h_GB);
E_0=0.424*sum(P_e_GT+P_e_pv+P_h_GB);
C4=0.75*(E_co2-E_0); %系统的碳交易成本
%总运行费用
C3=0.01*sum(P_e_GT+P_e_pv+P_h_GB); 
%% 目标函数
Obj=sum(pri_e.*P_buy)+3.5*sum(Gas)+0.01*sum(abs(P_e_tran))-0.03*sum(P_e_cut)+C3+C4...
    +0.016*sum(P_h_DR)+0.01*sum(abs(P_e_31)+abs(P_e_32))...
    +0.5*rou*(norm(P_e_31-z_31)^2)+sum(lambda_e_31.*P_e_31)...
    +0.5*rou*(norm(P_e_32-z_32)^2)+sum(lambda_e_32.*P_e_32);
%% 求解器配置与求解
ops=sdpsettings('solver','cplex','verbose',0,'usex0',0);
ops.cplex.mip.tolerances.mipgap=1e-6;
result=solvesdp(C,Obj,ops);
if result.problem == 0 
    P_e_31 = double(P_e_31);
    P_e_32 = double(P_e_32);
    Obj_MG3 = double(Obj);
    p3(1,:) = double(P_e_GT);
    p3(2,:) = double(P_e_pv);
    p3(3,:) = double(P_buy);
    p3(4,:) = double(P_batd);
    p3(5,:) = double(-1*P_batc);
    p3(6,:) = double(-1*L_e);
    p3(7,:) = -1*P_e_31;
    p3(8,:) = -1*P_e_32;
    p3(9,:) = double(-1*P_sell);
    h3(1,:) = double(P_h_GT);
    h3(2,:) = double(P_h_GB);
    h3(3,:) = double(-1*L_h);
    c3(1,1) = double(sum(pri_e.*P_buy));
    c3(1,2) = double(3.5*sum(Gas));
    c3(1,3) = double(sum(grid_sw.*P_sell));
    c3(1,4) = c3(1,1)+c3(1,2)-c3(1,3);
    c3(1,5) = double(0.01*sum(abs(P_e_tran))-0.03*sum(P_e_cut)+0.016*sum(P_h_DR));
    c3(1,6) = double(0.01*sum(abs(P_e_31)+abs(P_e_32)));
else
 disp('求解失败，失败原因为：')
 disp(result.info)
end
