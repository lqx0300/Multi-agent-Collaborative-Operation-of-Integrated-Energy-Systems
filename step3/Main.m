% 计及电能共享的基于非对称纳什谈判的多微网运行优化策略
%非对称纳什谈判,合作博弈,能量共济,多微网运行
%子问题2:微网间的非对称支付效益最大化

clc
clear
close all
U_no=[58217,65375,33627];
C_co=[40487,62688,31118];
load("Pe.mat");
P_e_12=Pe(1,:);
P_e_13=Pe(2,:);
P_e_23=Pe(3,:);
P_e_21=-P_e_12;
P_e_31=-P_e_13;
P_e_32=-P_e_23;
C_gap=U_no-C_co;
price_max = 100;
price_min=1;
pri_e=[0.40*ones(1,7),0.75*ones(1,4),1.20*ones(1,3),0.75*ones(1,4),1.20*ones(1,4),0.40*ones(1,2)];
%% ADMM迭代参数设置
%拉格朗日乘子初始化
lambda_e_12=zeros(1,24);%MG1和MG2之间的拉格朗日乘子
lambda_e_13=zeros(1,24);%MG1和MG3之间的拉格朗日乘子
lambda_e_21=zeros(1,24);%MG2和MG1之间的拉格朗日乘子
lambda_e_23=zeros(1,24);%MG2和MG3之间的拉格朗日乘子
lambda_e_31=zeros(1,24);%MG3和MG1之间的拉格朗日乘子
lambda_e_32=zeros(1,24);%MG3和MG2之间的拉格朗日乘子
maxIter=500;  %最大迭代次数
tolerant=1e-3;%收敛精度
rou = 1e-4;
iter=1;%迭代次数初始化
%对微网之间的交易量记录矩阵初始化
pri_e_12=zeros(maxIter+1,24);pri_e_21=zeros(maxIter+1,24);
pri_e_13=zeros(maxIter+1,24);pri_e_31=zeros(maxIter+1,24);
pri_e_23=zeros(maxIter+1,24);pri_e_32=zeros(maxIter+1,24);
kesai_12 = zeros(maxIter+1,24);
kesai_13 = zeros(maxIter+1,24);
kesai_21 = zeros(maxIter+1,24);
kesai_23 = zeros(maxIter+1,24);
kesai_31 = zeros(maxIter+1,24);
kesai_32 = zeros(maxIter+1,24);
r_list1 = [];
%% 迭代
while 1
    if iter==maxIter  %限制迭代次数
       disp('迭代不收敛,参数有误');
       break; 
    end 
    display(['迭代还未收敛,当前迭代第 ', num2str(iter),' 次']);
    [pri_e_12(iter+1,:),pri_e_13(iter+1,:),Obj_MG1(iter)]=Fun_MG1(kesai_12(iter,:),kesai_13(iter,:),lambda_e_12,lambda_e_13,rou,price_max,price_min,P_e_12,P_e_13,C_gap(1));
    [pri_e_21(iter+1,:),pri_e_23(iter+1,:),Obj_MG2(iter)]=Fun_MG2(kesai_21(iter,:),kesai_23(iter,:),lambda_e_21,lambda_e_23,rou,price_max,price_min,P_e_21,P_e_23,C_gap(2));
    [pri_e_31(iter+1,:),pri_e_32(iter+1,:),Obj_MG3(iter)]=Fun_MG3(kesai_31(iter,:),kesai_32(iter,:),lambda_e_31,lambda_e_32,rou,price_max,price_min,P_e_31,P_e_32,C_gap(3));
    kesai_12(iter+1,:) = 0.5*(pri_e_12(iter+1,:)+pri_e_21(iter+1,:));
    kesai_13(iter+1,:) = 0.5*(pri_e_13(iter+1,:)+pri_e_31(iter+1,:));
    kesai_21(iter+1,:) = 0.5*(pri_e_21(iter+1,:)+pri_e_12(iter+1,:));
    kesai_23(iter+1,:) = 0.5*(pri_e_23(iter+1,:)+pri_e_32(iter+1,:));
    kesai_31(iter+1,:) = 0.5*(pri_e_31(iter+1,:)+pri_e_13(iter+1,:));
    kesai_32(iter+1,:) = 0.5*(pri_e_32(iter+1,:)+pri_e_23(iter+1,:));
    lambda_e_12=lambda_e_12+rou*(pri_e_12(iter+1,:)-kesai_12(iter+1,:));
    lambda_e_13=lambda_e_13+rou*(pri_e_13(iter+1,:)-kesai_13(iter+1,:));
    lambda_e_21=lambda_e_21+rou*(pri_e_21(iter+1,:)-kesai_21(iter+1,:));
    lambda_e_23=lambda_e_23+rou*(pri_e_23(iter+1,:)-kesai_23(iter+1,:));
    lambda_e_31=lambda_e_31+rou*(pri_e_31(iter+1,:)-kesai_31(iter+1,:));
    lambda_e_32=lambda_e_32+rou*(pri_e_32(iter+1,:)-kesai_32(iter+1,:));
    r = sum(norm(pri_e_12(iter+1,:)-kesai_12(iter+1,:))^2)+sum(norm(pri_e_13(iter+1,:)-kesai_13(iter+1,:))^2);
    s = rou*(sum(norm(kesai_12(iter+1,:)-kesai_12(iter,:))^2)+sum(norm(kesai_13(iter+1,:)-kesai_13(iter,:))^2));
    r_list1 = [r_list1,r];
    if sqrt(r)>10*sqrt(s)
        rou = 2*rou;
    elseif sqrt(s)>10*sqrt(r)
        rou = 0.5*rou;
    end
    %判断收敛条件
    if r<=tolerant
       display(['迭代收敛,在第 ', num2str(iter),' 次收敛']);
       break; 
    end
    iter=iter+1;
end

