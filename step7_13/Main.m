%% 计及电能共享的基于非对称纳什谈判的多微网运行优化策略
%非对称纳什谈判,合作博弈,能量共济,多微网运行
%子问题1:微网间的社会成本最小化问题


clc
clear
close all

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
iter=1;%迭代次数初始化
rou = 1e-4;
U_no=[58217,65375,33627];
%对微网之间的交易量记录矩阵初始化
P_e_12=zeros(maxIter+1,24);P_e_21=zeros(maxIter+1,24);
P_e_13=zeros(maxIter+1,24);P_e_31=zeros(maxIter+1,24);
P_e_23=zeros(maxIter+1,24);P_e_32=zeros(maxIter+1,24);
z_12 = zeros(maxIter+1,24);z_21 = zeros(maxIter+1,24);
z_13 = zeros(maxIter+1,24);z_31 = zeros(maxIter+1,24);
z_23 = zeros(maxIter+1,24);z_32 = zeros(maxIter+1,24);
%记录

%% 迭代
while 1
    if iter==maxIter  %限制迭代次数
       disp('迭代不收敛,参数有误');
       break; 
    end 
    display(['迭代还未收敛,当前迭代第 ', num2str(iter),' 次']);
    [P_e_12(iter+1,:),P_e_13(iter+1,:),Obj_MG1(iter),p1,h1,c1]=Fun_MG1(z_12(iter,:),z_13(iter,:),lambda_e_12,lambda_e_13,rou);
    [P_e_21(iter+1,:),P_e_23(iter+1,:),Obj_MG2(iter),p2,h2,c2]=Fun_MG2(z_21(iter,:),z_23(iter,:),lambda_e_21,lambda_e_23,rou);
    [P_e_31(iter+1,:),P_e_32(iter+1,:),Obj_MG3(iter),p3,h3,c3]=Fun_MG3(z_31(iter,:),z_32(iter,:),lambda_e_31,lambda_e_32,rou);
    z_12(iter+1,:)=0.5*(P_e_12(iter+1,:)-P_e_21(iter+1,:));
    z_13(iter+1,:)=0.5*(P_e_13(iter+1,:)-P_e_31(iter+1,:));
    z_21(iter+1,:)=0.5*(P_e_21(iter+1,:)-P_e_12(iter+1,:));
    z_23(iter+1,:)=0.5*(P_e_23(iter+1,:)-P_e_32(iter+1,:));
    z_31(iter+1,:)=0.5*(P_e_31(iter+1,:)-P_e_13(iter+1,:));
    z_32(iter+1,:)=0.5*(P_e_32(iter+1,:)-P_e_23(iter+1,:));
    lambda_e_12=lambda_e_12+rou*(P_e_12(iter+1,:)-z_12(iter+1,:));
    lambda_e_13=lambda_e_13+rou*(P_e_13(iter+1,:)-z_13(iter+1,:));
    lambda_e_21=lambda_e_21+rou*(P_e_21(iter+1,:)-z_21(iter+1,:));
    lambda_e_23=lambda_e_23+rou*(P_e_23(iter+1,:)-z_23(iter+1,:));
    lambda_e_31=lambda_e_31+rou*(P_e_31(iter+1,:)-z_31(iter+1,:));
    lambda_e_32=lambda_e_32+rou*(P_e_32(iter+1,:)-z_32(iter+1,:));
    r = sum(norm(P_e_12(iter+1,:)-z_12(iter+1,:))^2)+sum(norm(P_e_13(iter+1,:)-z_13(iter+1,:))^2);
    % s = rou*(sum(norm(z_12(iter+1,:)-z_12(iter,:))^2)+sum(norm(z_13(iter+1,:)-z_13(iter,:))^2));
    % if sqrt(r)>10*sqrt(s)
    %     rou = 2*rou;
    % elseif sqrt(s)>10*sqrt(r)
    %     rou = 0.5*rou;
    % end
    %判断收敛条件
    if (r<=tolerant)|(iter==200) 
       display(['迭代收敛,在第 ', num2str(iter),' 次收敛']);
       break; 
    end
    iter=iter+1;
end
