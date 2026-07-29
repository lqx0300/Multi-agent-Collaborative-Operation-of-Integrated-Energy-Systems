function [P_e_21,P_e_23,Obj_MG2,p2,h2,c2]=Fun_MG2(z_21,z_23,lambda_e_21,lambda_e_23,rou)
%微网2(MG2)的分布式优化迭代模型
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
P_h_GB=sdpvar(1,24); %余热锅炉的产热功率
P_buy=sdpvar(1,24);  %微网向外电网的购买的电功率
P_sell=sdpvar(1,24); %微网向外电网的售出的电功率
Gas_GT=sdpvar(1,24); %GT的耗气量
Gas_GB=sdpvar(1,24); %GB的耗气量
Gas=sdpvar(1,24);    %系统的总耗气量
%P2G+CCS
P_e1=sdpvar(1,24); %CHP的供电功率
P_e2=sdpvar(1,24); %CHP的供给P2G的功率
P_e3=sdpvar(1,24); %CHP的供给CCS的功率
P_h=sdpvar(1,24);  %CHP的输出热功率
P_gs=sdpvar(1,24); %P2G的产气功率
C_cc=sdpvar(1,24); %CCS的碳捕集量/P2G所用的二氧化碳量
P_e_21=sdpvar(1,24);
P_e_23=sdpvar(1,24);
p2=zeros(9,24);
h2=zeros(3,24);
c2=zeros(1,7);
%% 导入电/热负荷和电网购电电价
L_e0=[1774,1450,1296,1219,1095,1265,1481,1944,2484,2083,1651,1188,1080,1126,1033,1033,941,1450,2283,3148,3904,3719,2746,2453]; 
L_h0=[1610,1594,1594,1610,1633,1633,1286,1201,1117,1109,1648,1656,1664,1140,1124,1109,1286,1309,1302,1325,1479,1502,1340,1332]*0.5;
Predict_pv=[0,0,0,0,0,0,967,1287,1583,1833,1918,1942,2004,1957,1669,1076,655,0,0,0,0,0,0,0]*2;
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
%带P2G和CCS的CHP运行约束
C=[C,
    0-P_e2-P_e3<=P_e1<=3000-P_e2-P_e3, %CHP的供电功率约束
    0<=P_e2<=300, %P2G设备的耗电功率约束
    0<=P_e3<=600, %CCS设备的耗电功率约束
    0<=P_e1<=2000, %CHP的供电功率上下限约束,公式(11)
    0<=P_e1, %CHP的供电功率非负性约束
    max((0-0.15*P_h-P_e2-P_e3),(0.85*(P_h)-P_e2-P_e3))<=P_e1<=3000-0.20*P_h-P_e2-P_e3, %CHP的热电耦合约束    式（13）
    max((0-0.15*P_h),(0.85*(P_h-50)-300-600))<=P_e1<=3000-0.20*P_h-0-0, %考虑P2G和CCS后的CHP的热电耦合约束       式（15）
    (0.55/(1+0.5*1.02))*max((0-0.15*P_h-P_e1),(0.85*(P_h-50)-P_e1))<=P_gs<=(0.55/(1+0.5*1.02))*(3000-0.20*P_h-P_e1), %产气功率上下限约束  式（17）
    -800<=(P_e1(2:24)+P_e2(2:24)+P_e3(2:24))-(P_e1(1:23)+P_e2(1:23)+P_e3(1:23))<=800, %CHP的爬坡约束   
    P_gs==0.55*P_e2, %P2G产气功率与耗电量约束  式（2）
    C_cc==1.02*P_e2, %P2G运行所需要的二氧化碳量与电功率约束  式（3）
    P_e3==0.55*C_cc, %CCS的耗电量与碳捕集量约束   式（4）
  ];
%CCS的最大碳捕集量
C=[C,0<=C_cc<=0.55*(P_e1+P_e2+P_e3+0.15*P_h),]; %此式为全体的碳排放量，最大捕集量不会超过最大排放量
%CHP机组以及GB等设备运行约束
for t=1:24
    C=[C,
       0<=P_h_GB(t)<=500, %余热锅炉产热功率上下限约束
       0<=P_e_pv(t)<=Predict_pv(t), %风力发电上下限约束
      ];
end
%微网的热/电负荷平衡约束
for t=1:24
    C=[C,
       P_e1(t)+P_e_pv(t)+P_buy(t)+P_batd(t)==P_batc(t)+L_e(t)+P_e_21(t)+P_e_23(t)+P_sell(t),
       P_h(t)+P_h_GB(t)==L_h(t),
      ];
end
%变量非负性等约束
for t=1:24
    C=[C,
       P_buy(t)>=0,
       0<=P_sell(t)<=2000,
       2000>=P_e_21(t)>=-2000,
       2000>=P_e_23(t)>=-2000,
      ];
end
%
for t=1:24
    C=[C,
       P_e1(t)+P_e2(t)+P_e3(t)==0.35*Gas_GT(t), %GT耗气量约束
       P_h_GB(t)==0.9*Gas_GB(t), %GB耗气量约束
       Gas(t)==Gas_GT(t)+Gas_GB(t)-P_gs(t),%总耗气量约束
      ];
end
%碳交易部分
E_co2=0.55*sum(P_e1+P_e2+P_e3+P_h)+0.65*sum(P_h_GB)-sum(C_cc);
E_0=0.424*sum(P_e1+P_e2+P_e3+P_e_pv+P_h_GB);
C4=0.75*(E_co2-E_0); %系统的碳交易成本
%总运行费用
C3=0.01*sum(P_e1+P_e2+P_e3+P_e_pv+P_h_GB);
%% 目标函数
Obj=sum(pri_e.*P_buy)+3.5*sum(Gas)+0.01*sum(abs(P_e_tran))-0.03*sum(P_e_cut)+C3+C4...
    +0.016*sum(P_h_DR)+0.01*sum(abs(P_e_21)+abs(P_e_23))...
    +0.5*rou*(norm(P_e_21-z_21)^2)+sum(lambda_e_21.*P_e_21)...
    +0.5*rou*(norm(P_e_23-z_23)^2)+sum(lambda_e_23.*P_e_23);
%% 求解器配置与求解
ops=sdpsettings('solver','cplex','verbose',0,'usex0',0);
ops.cplex.mip.tolerances.mipgap=1e-6;
result=solvesdp(C,Obj,ops);
if result.problem == 0 
    P_e_21 = double(P_e_21);
    P_e_23 = double(P_e_23);
    Obj_MG2 = double(Obj);
    p2(1,:) = double(P_e1);
    p2(2,:) = double(P_e_pv);
    p2(3,:) = double(P_buy);
    p2(4,:) = double(P_batd);
    p2(5,:) = double(-1*P_batc);
    p2(6,:) = double(-1*L_e);
    p2(7,:) = -1*P_e_21;
    p2(8,:) = -1*P_e_23;
    p2(9,:) = double(-1*P_sell);
    h2(1,:) = double(P_h);
    h2(2,:) = double(P_h_GB);
    h2(3,:) = double(-1*L_h);
    c2(1,1) = double(sum(pri_e.*P_buy));
    c2(1,2) = double(3.5*sum(Gas));
    c2(1,3) = double(sum(grid_sw.*P_sell));
    c2(1,4) = c2(1,1)+c2(1,2)-c2(1,3);
    c2(1,5) = double(0.01*sum(abs(P_e_tran))-0.03*sum(P_e_cut)+0.016*sum(P_h_DR));
    c2(1,6) = double(0.01*sum(abs(P_e_21)+abs(P_e_23)));
    c2(1,7) = C4;
else
 disp('求解失败，失败原因为：')
 disp(result.info)
end
