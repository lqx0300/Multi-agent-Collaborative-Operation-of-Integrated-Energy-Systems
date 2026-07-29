clc;
clear;
U_no=[58217,65375,33627];
C_co=[40487,62688,31118];
c_gap = [20993,3396,583];
kesai = zeros(100,3);
kesai_max = zeros(1,3);
c_gap_lie = zeros(32,3);
c_cap_sum = zeros(1,32);
eps = 0.0001;
sigma = [0.1,0.1,0.1];
for i = 1:3
    kesai(1,i) = 0.0001;
end
iter = 1;
maxIter = 100;
while 1
    if iter==maxIter  %限制迭代次数
       disp('迭代不收敛,参数有误');
       break; 
    end 
    display(['迭代还未收敛,当前迭代第 ', num2str(iter),' 次']);
    for i = 1:3
        c_gap_lie(iter,i) = c_gap(i)-kesai(iter,i)*abs(c_gap(i));
    end
    c_cap_sum(iter) = sum(c_gap_lie(iter,:));
    for i = 1:3
        kesai_max(i) = (c_gap(i) + sum(c_gap_lie(iter,:))-c_gap_lie(iter,i))/abs(c_gap(i));
    end
    % for i = 1:3
    %     sigma(i) = abs((kesai_max(i)-kesai(i))/(4*kesai_max(i)));
    % end
    for i = 1:3
        kesai(iter+1,i) = (1-sigma(i))*kesai(iter,i)+sigma(i)*kesai_max(i);
    end
    if sum(kesai(iter+1,:)-kesai(iter,:)) <= eps
       display(['迭代收敛,在第 ', num2str(iter),' 次收敛']);
       break; 
    end
    iter=iter+1;
end
c_gap_lie_1 = c_gap(1) - kesai(iter+1,1)*abs(c_gap(1))
c_gap_lie_2 = c_gap(2) - kesai(iter+1,2)*abs(c_gap(2))
c_gap_lie_3 = c_gap(3) - kesai(iter+1,3)*abs(c_gap(3))
c_gap_list = zeros(1,32)';
c_gap_list = [c_gap_list,c_gap_lie,c_cap_sum'];
save("c_gap","c_gap_list")




