%%%%%%%%%%%%%%%%%%%%%%%%%%%%%用户分组第二阶段%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nums_cell1=0;
for i=1:nums_ue
    if(UE_information(i,2)==1)
        nums_cell1=nums_cell1+1;
    end
end
cell1=cell(nums_cell1,9);%第九列放W
j=0;
for i=1:nums_ue
    if(UE_information(i,2)==1)
        j=j+1;
        numsrows=size(W{i,1},1);
        cell1{j,9}=W{i,1}(1:numsrows/2, :);%%%%%
        cell1{j,1}=j;
        cell1{j,2}=H{i,1};%第二列是信道矩阵
        cell1{j,7}=theta2(i,1);%%第七列放一下小区间干扰
    end
end
for i=1:nums_cell1
    cell1{i,3}=(cell1{i,2}')*cell1{i,2};%第三列所示H'H
end
for i=1:nums_cell1
    [V,D]=eig(cell1{i,3});
    %disp('特征值矩阵 D:');
    %disp(D);
    cell1{i,4}=D;
    %第四列就是lamda

    %disp('特征向量矩阵 V:');
    %disp(V);
    cell1{i,5}=V;
%第五列是xi
end
for i=1:nums_cell1
 tao=zeros(Mt,1);
    for j=1:Mt
        tao=cell1{i,4}(j,j)*cell1{i,5}(j,:);
    end
    cell1{i,6}=tao;
end

row=zeros(nums_cell1,nums_cell1);%信道相关性矩阵
for i=1:nums_cell1
    for j=1:nums_cell1
        if i==j
            row(i,j)=0;
        end
        row(i,j)=norm(((cell1{i,6})'*cell1{i,6})'*((cell1{j,6})'*cell1{j,6}))/(norm(((cell1{i,6})'*cell1{i,6}))*norm(((cell1{j,6})'*cell1{j,6}))); 
    end
end
v=pi^2*cos(0)*(Mr-1)/12;
B=100000000;%带宽总和
P=10;%每个基站传输功率
n0=-174;
n=10^(n0/10)/1000;

aa=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
result_1=[];a=0;
result_MI_1=[];
for a_x=0:10
      a=0.1*a_x;
  %a=0.5;

nums_group=0;
group_information=cell(100,100);%第一列是组数，第二列是成员信息
%将第一个用户加入组1
nums_group=nums_group+1;
group_information{1,1}=nums_group;
group1=[];group1(1,1)=1;
group1_nums=1;
group_matrix=zeros(nums_cell1,1);%用来标记用户是否进行了分组
group_matrix(1,1)=1;
for j=2:nums_cell1
    if row(1,j)<a
        group1_nums=group1_nums+1;
        group1(group1_nums,1)=j;
        group_matrix(j,1)=1;
    end
end
group_information{nums_group,2}=group1;
%%第二阶段步骤3
%for jishu=1:nums_cell1
    while(any(group_matrix==0))
        countzeros=sum(group_matrix==0);
        nums_group=nums_group+1;
        [startindex,colindex]=find(group_matrix==0);
        group_new=[]; 
        group_new(1,1)=startindex(1,1);
        nums_group_new=1;
        group_matrix(startindex(1))=1;
        j=countzeros;
        [startindex_find,colindex_find]=find(group_matrix==0);%找到另外的还没分组的用户
        for i=1:j-1
            if row(startindex(1,1),startindex_find(i,1))<a
               nums_group_new=nums_group_new+1;
               group_new(nums_group_new,1)=startindex_find(i);
               group_matrix(startindex_find(i,1))=1;
            end
        end
    group_information{nums_group,2}=group_new;
    j=0;
    clear group_new;
    end
   
%end
for i=1:nums_group
    group_information{i,1}=i;
end


%%第二阶段步骤4
for i=1:nums_group
    nums_group_i=length(group_information{i,2});
    for j=1:nums_group_i
        nums_less_than_a=0;
        for k=j+1:nums_group_i
            if k==j
                continue
            end
            if row(group_information{i,2}(j,1),group_information{i,2}(k,1))>a
                nums_less_than_a=nums_less_than_a+1;
            end
        end
        if nums_less_than_a>1
            %添加新组
            nums_group=nums_group+1;
            group_new=[];
            group_new=group_information{i,2}(j);
            group_information{nums_group,1}=nums_group;
            group_information{nums_group,2}=group_new;
            nums_group_i= nums_group_i-1;%%回溯防止越界
            %从该分组中剔除该用户
            group_information{i,2}(j)=[];

        end
    end
end
%disp(nums_group);
%%%%%%%%%%%%SINR计算%%%%%%%%%%
%%此阶段不考虑资源分配的问题
%%%每个小区 每个用户平均分配功率， P=5w
%%%每个小区 总带宽B,
%SINR=P/theta2+噪声
P_every_ue=P/nums_cell1;
every_group=cell(nums_group,3);%第一列组内用户数量，第二列组的带宽,第三列用户成员
for i=1:nums_group
    every_group{i,1}=length(group_information{i,2});
    every_group{i,2}=B*(length(group_information{i,2})/nums_cell1);%组内复用带宽，组间带宽正交
    every_group{i,3}=group_information{i,2};
end
R_sum=0;
MI_sum=0;
H_cell1=[];
trace1=0;
for i=1:nums_cell1
    H_cell1=[H_cell1,cell1{i,2}'];
end
H_cell1=(H_cell1')*H_cell1;
for i=1:length(H_cell1(1,:))
    for j=1:length(H_cell1(:,1))
        if i==j
        trace1=trace1+H_cell1(i,j);
        end
    end
end

some111=zeros(nums_group,100);

for i=1:nums_group%计算R
    for j=1:every_group{i,1}
        hk=cell1{every_group{i,3}(j,1),2};
        wk=cell1{every_group{i,3}(j,1),9};
        disturb=0;
        for l=1:every_group{i,1}
            if l==j
                continue
            end
            hl=cell1{every_group{i,3}(l,1),2};
            wl=cell1{every_group{i,3}(l,1),9};
            disturb=disturb+norm(wl'*hk*hl'/sqrt(trace1))*P_every_ue;
        end
  
        %%sssssss=P_every_ue*guiyihuaxishu/sqrt(trace1);
       
        fenzi=norm(wk'*(hk)*hk'/sqrt(trace1));
        SINR=P_every_ue*fenzi/( disturb+every_group{i,2}*n+2*P*cell1{every_group{i,3}(j,1),7});
        some=log2(1+SINR);
        %some111(i,j)=some;
        R=every_group{i,2}*some;
        R_sum=R_sum+R;
    end
end
%disp(nums_group);
%disp(R_sum);
 result_1(a_x+1,1)=R_sum;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%用户分组第二阶段%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nums_cell1=0;
for i=1:nums_ue
    if(UE_information(i,2)==2)
        nums_cell1=nums_cell1+1;
    end
end
cell1=cell(nums_cell1,9);%第九列放W
j=0;
for i=1:nums_ue
    if(UE_information(i,2)==2)
        j=j+1;
        numsrows=size(W{i,1},1);
        cell1{j,9}=W{i,1}(1:numsrows/2, :);%%%%%
        cell1{j,1}=j;
        cell1{j,2}=H{i,1};%第二列是信道矩阵
        cell1{j,7}=theta2(i,1);%%第七列放一下小区间干扰
    end
end
for i=1:nums_cell1
    cell1{i,3}=(cell1{i,2}')*cell1{i,2};%第三列所示H'H
end
for i=1:nums_cell1
    [V,D]=eig(cell1{i,3});
    %disp('特征值矩阵 D:');
    %disp(D);
    cell1{i,4}=D;
    %第四列就是lamda

    %disp('特征向量矩阵 V:');
    %disp(V);
    cell1{i,5}=V;
%第五列是xi
end
for i=1:nums_cell1
 tao=zeros(Mt,1);
    for j=1:Mt
        tao=cell1{i,4}(j,j)*cell1{i,5}(j,:);
    end
    cell1{i,6}=tao;
end

row=zeros(nums_cell1,nums_cell1);%信道相关性矩阵
for i=1:nums_cell1
    for j=1:nums_cell1
        if i==j
            row(i,j)=0;
        end
        row(i,j)=norm(((cell1{i,6})'*cell1{i,6})'*((cell1{j,6})'*cell1{j,6}))/(norm(((cell1{i,6})'*cell1{i,6}))*norm(((cell1{j,6})'*cell1{j,6}))); 
    end
end
v=pi^2*cos(0)*(Mr-1)/12;
B=100000000;%带宽总和
P=10;%每个基站传输功率
n0=-174;
n=10^(n0/10)/1000;
disturb=0.1;
D=40;
omega=2*pi;

aa=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
result_2=[];a=0;
result_MI_2=[];
for a_x=0:10
      a=0.1*a_x;%test 分组阈值

  %a=0.5;

nums_group=0;
group_information=cell(100,100);%第一列是组数，第二列是成员信息
%将第一个用户加入组1
nums_group=nums_group+1;
group_information{1,1}=nums_group;
group1=[];group1(1,1)=1;
group1_nums=1;
group_matrix=zeros(nums_cell1,1);%用来标记用户是否进行了分组
group_matrix(1,1)=1;
for j=2:nums_cell1
    if row(1,j)<a
        group1_nums=group1_nums+1;
        group1(group1_nums,1)=j;
        group_matrix(j,1)=1;
    end
end
group_information{nums_group,2}=group1;
%%第二阶段步骤3
%for jishu=1:nums_cell1
    while(any(group_matrix==0))
        countzeros=sum(group_matrix==0);
        nums_group=nums_group+1;
        [startindex,colindex]=find(group_matrix==0);
        group_new=[]; 
        group_new(1,1)=startindex(1,1);
        nums_group_new=1;
        group_matrix(startindex(1))=1;
        j=countzeros;
        [startindex_find,colindex_find]=find(group_matrix==0);%找到另外的还没分组的用户
        for i=1:j-1
            if row(startindex(1,1),startindex_find(i,1))<a
               nums_group_new=nums_group_new+1;
               group_new(nums_group_new,1)=startindex_find(i);
               group_matrix(startindex_find(i,1))=1;
            end
        end
    group_information{nums_group,2}=group_new;
    j=0;
    clear group_new;
    end
   
%end
for i=1:nums_group
    group_information{i,1}=i;
end

%%第二阶段步骤4
for i=1:nums_group
    nums_group_i=length(group_information{i,2});
    for j=1:nums_group_i
        nums_less_than_a=0;
        for k=j+1:nums_group_i
            if k==j
                continue
            end
            if row(group_information{i,2}(j,1),group_information{i,2}(k,1))>a
                nums_less_than_a=nums_less_than_a+1;
            end
        end
        if nums_less_than_a>1
            %添加新组
            nums_group=nums_group+1;
            group_new=[];
            group_new=group_information{i,2}(j);
            group_information{nums_group,1}=nums_group;
            group_information{nums_group,2}=group_new;
            nums_group_i= nums_group_i-1;%%回溯防止越界
            %从该分组中剔除该用户
            group_information{i,2}(j)=[];

        end
    end
end
%disp(nums_group);
%%%%%%%%%%%%SINR计算%%%%%%%%%%
%%此阶段不考虑资源分配的问题
%%%每个小区 每个用户平均分配功率， P=5w
%%%每个小区 总带宽B,
%SINR=P/theta2+噪声
P_every_ue=P/nums_cell1;
every_group=cell(nums_group,3);%第一列组内用户数量，第二列组的带宽,第三列用户成员
for i=1:nums_group
    every_group{i,1}=length(group_information{i,2});
    every_group{i,2}=B*(length(group_information{i,2})/nums_cell1);%组内复用带宽，组间带宽正交
    every_group{i,3}=group_information{i,2};
end
R_sum=0;
MI_sum=0;
H_cell1=[];
trace1=0;
for i=1:nums_cell1
    H_cell1=[H_cell1,cell1{i,2}'];
end
H_cell1=(H_cell1')*H_cell1;
for i=1:length(H_cell1(1,:))
    for j=1:length(H_cell1(:,1))
        if i==j
        trace1=trace1+H_cell1(i,j);
        end
    end
end

some111=zeros(nums_group,100);

for i=1:nums_group%计算R
    for j=1:every_group{i,1}
        hk=cell1{every_group{i,3}(j,1),2};
        wk=cell1{every_group{i,3}(j,1),9};
        disturb=0;
        for l=1:every_group{i,1}
            if l==j
                continue
            end
            hl=cell1{every_group{i,3}(l,1),2};
            wl=cell1{every_group{i,3}(l,1),9};
            disturb=disturb+norm(wl'*hk*hl'/sqrt(trace1))*P_every_ue;
        end
  
        %%sssssss=P_every_ue*guiyihuaxishu/sqrt(trace1);
        fenzi=norm(wk'*(hk)*hk'/sqrt(trace1));
        SINR=P_every_ue*fenzi/( disturb+every_group{i,2}*n+2*P*cell1{every_group{i,3}(j,1),7});
        some=log2(1+SINR);
        %some111(i,j)=some;
        R=every_group{i,2}*some;
        R_sum=R_sum+R;
    end
end
%disp(nums_group);
%disp(R_sum);
 result_2(a_x+1,1)=R_sum;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%小区三%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nums_cell1=0;
for i=1:nums_ue
    if(UE_information(i,2)==3)
        nums_cell1=nums_cell1+1;
    end
end
cell1=cell(nums_cell1,9);%第九列放W
j=0;
for i=1:nums_ue
    if(UE_information(i,2)==3)
        j=j+1;
        numsrows=size(W{i,1},1);
        cell1{j,9}=W{i,1}(1:numsrows/2, :);%%%%%
        cell1{j,1}=j;
        cell1{j,2}=H{i,3};%第二列是信道矩阵
        cell1{j,7}=theta2(i,3);%%第七列放一下小区间干扰
    end
end
for i=1:nums_cell1
    cell1{i,3}=(cell1{i,2}')*cell1{i,2};%第三列所示H'H
end
for i=1:nums_cell1
    [V,D]=eig(cell1{i,3});
    %disp('特征值矩阵 D:');
    %disp(D);
    cell1{i,4}=D;
    %第四列就是lamda

    %disp('特征向量矩阵 V:');
    %disp(V);
    cell1{i,5}=V;
%第五列是xi
end
for i=1:nums_cell1
 tao=zeros(Mt,1);
    for j=1:Mt
        tao=cell1{i,4}(j,j)*cell1{i,5}(j,:);
    end
    cell1{i,6}=tao;
end

row=zeros(nums_cell1,nums_cell1);%信道相关性矩阵
for i=1:nums_cell1
    for j=1:nums_cell1
        if i==j
            row(i,j)=0;
        end
        row(i,j)=norm(((cell1{i,6})'*cell1{i,6})'*((cell1{j,6})'*cell1{j,6}))/(norm(((cell1{i,6})'*cell1{i,6}))*norm(((cell1{j,6})'*cell1{j,6}))); 
    end
end
B=100000000;%带宽总和
P=10;%每个基站传输功率
n0=-174;
n=10^(n0/10)/1000;
disturb=0.1;

aa=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
result_3=[];a=0;
result_MI_3=[];
for a_x=0:10
      a=0.1*a_x;%test 分组阈值
  %a=0.5;
nums_group=0;
group_information=cell(100,100);%第一列是组数，第二列是成员信息
%将第一个用户加入组1
nums_group=nums_group+1;
group_information{1,1}=nums_group;
group1=[];group1(1,1)=1;
group1_nums=1;
group_matrix=zeros(nums_cell1,1);%用来标记用户是否进行了分组
group_matrix(1,1)=1;
for j=2:nums_cell1
    if row(1,j)<a
        group1_nums=group1_nums+1;
        group1(group1_nums,1)=j;
        group_matrix(j,1)=1;
    end
end
group_information{nums_group,2}=group1;
%%第二阶段步骤3
%for jishu=1:nums_cell1
    while(any(group_matrix==0))
        countzeros=sum(group_matrix==0);
        nums_group=nums_group+1;
        [startindex,colindex]=find(group_matrix==0);
        group_new=[]; 
        group_new(1,1)=startindex(1,1);
        nums_group_new=1;
        group_matrix(startindex(1))=1;
        j=countzeros;
        [startindex_find,colindex_find]=find(group_matrix==0);%找到另外的还没分组的用户
        for i=1:j-1
            if row(startindex(1,1),startindex_find(i,1))<a
               nums_group_new=nums_group_new+1;
               group_new(nums_group_new,1)=startindex_find(i);
               group_matrix(startindex_find(i,1))=1;
            end
        end
    group_information{nums_group,2}=group_new;
    j=0;
    clear group_new;
    end
   
%end
for i=1:nums_group
    group_information{i,1}=i;
end

%%第二阶段步骤4
for i=1:nums_group
    nums_group_i=length(group_information{i,2});
    for j=1:nums_group_i
        nums_less_than_a=0;
        for k=j+1:nums_group_i
            if k==j
                continue
            end
            if row(group_information{i,2}(j,1),group_information{i,2}(k,1))>a
                nums_less_than_a=nums_less_than_a+1;
            end
        end
        if nums_less_than_a>1
            %添加新组
            nums_group=nums_group+1;
            group_new=[];
            group_new=group_information{i,2}(j);
            group_information{nums_group,1}=nums_group;
            group_information{nums_group,2}=group_new;
            nums_group_i= nums_group_i-1;%%回溯防止越界
            %从该分组中剔除该用户
            group_information{i,2}(j)=[];
        end
    end
end
%disp(nums_group);
%%%%%%%%%%%%SINR计算%%%%%%%%%%
%%此阶段不考虑资源分配的问题
%%%每个小区 每个用户平均分配功率， P=5w
%%%每个小区 总带宽B,
%SINR=P/theta2+噪声
P_every_ue=P/nums_cell1;
every_group=cell(nums_group,3);%第一列组内用户数量，第二列组的带宽,第三列用户成员
for i=1:nums_group
    every_group{i,1}=length(group_information{i,2});
    every_group{i,2}=B*(length(group_information{i,2})/nums_cell1);%组内复用带宽，组间带宽正交
    every_group{i,3}=group_information{i,2};
end
R_sum=0;
MI_sum=0;
H_cell1=[];
trace1=0;
for i=1:nums_cell1
    H_cell1=[H_cell1,cell1{i,2}'];
end
H_cell1=(H_cell1')*H_cell1;
for i=1:length(H_cell1(1,:))
    for j=1:length(H_cell1(:,1))
        if i==j
        trace1=trace1+H_cell1(i,j);
        end
    end
end
some111=zeros(nums_group,100);
for i=1:nums_group%计算R
    for j=1:every_group{i,1}
        hk=cell1{every_group{i,3}(j,1),2};
        wk=cell1{every_group{i,3}(j,1),9};
        disturb=0;
        for l=1:every_group{i,1}
            if l==j
                continue
            end
            hl=cell1{every_group{i,3}(l,1),2};
            wl=cell1{every_group{i,3}(l,1),9};
            disturb=disturb+norm(wl'*hk*hl'/sqrt(trace1))*P_every_ue;
        end
        %%sssssss=P_every_ue*guiyihuaxishu/sqrt(trace1);
        fenzi=norm(wk'*(hk)*hk'/sqrt(trace1));
        SINR=P_every_ue*fenzi/( disturb+every_group{i,2}*n+2*P*cell1{every_group{i,3}(j,1),7});
        some=log2(1+SINR);
        %some111(i,j)=some;
        R=every_group{i,2}*some;
        R_sum=R_sum+R;
    end
end
%disp(nums_group);
%disp(R_sum);
 result_3(a_x+1,1)=R_sum;
end

figure;
   result_1=result_1';
 values1 = spcrv([[aa(1) aa aa(end)];[result_1(1) result_1 result_1(end)]],4);
plot(values1(1,:),values1(2,:), '-k','linewidth',2);

hold on;
 result_2=result_2';
 values2 = spcrv([[aa(1) aa aa(end)];[result_2(1) result_2 result_2(end)]],4);
plot(values2(1,:),values2(2,:), '-r','linewidth',2);
result_3=result_3';
 values3 = spcrv([[aa(1) aa aa(end)];[result_3(1) result_3 result_3(end)]],4);
 plot(values3(1,:),values3(2,:), '-b','linewidth',2);

plot(values1(1,2:16:end),values1(2,2:16:end), 'k o','linewidth',2);
  plot(values2(1,2:16:end),values2(2,2:16:end), '+ r ','linewidth',2);
  plot(values3(1,2:16:end),values3(2,2:16:end), '< b ','linewidth',2);
  legend('BS1','BS2','BS3');
xlim([0.25, 0.95])
grid on;
xlabel('grouping threshold \alpha');
ylabel('total data rate(bits/s)');