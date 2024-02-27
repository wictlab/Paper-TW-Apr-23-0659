clc;
clear all;

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

%%create a square which len is 120 meter
N = 30;  %radius
%%three base station
BS = zeros(3,2);
BS(1,1) = 10 + 15*sqrt(3);
BS(1,2) =  62;
BS(2,1) = BS(1,1) + 45;
BS(2,2) = BS(1,2) + 15*sqrt(3);
BS(3,1) = BS(2,1);
BS(3,2) = BS(1,2) - 15*sqrt(3);
plot(BS(1,1),BS(1,2),'*','color','b','markersize',10);
hold on;
plot(BS(2,1),BS(2,2),'*','color','m','markersize',10);
hold on;
plot(BS(3,1),BS(3,2),'*','color','r','markersize',10);
hold on;

%%plot the domain of each base station
A(1) = BS(1,1);A(2) = BS(1,2);
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'b','LineWidth',2,'HandleVisibility','off')
hold on;
A(1) = BS(2,1);A(2) = BS(2,2);
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'m','LineWidth',2,'HandleVisibility','off')
hold on;
A(1) = BS(3,1);A(2) = BS(3,2);
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'r','LineWidth',2,'HandleVisibility','off')
hold on;

%%plot users
%% 12 users in each cell
%first cell
u1 = zeros(10,2);
u1(1,1) = BS(1,1) - 20;
u1(1,2) = BS(1,2);
u1(2,1) = u1(1,1) + 2.5;
u1(2,2) = u1(1,2) + 4;
u1(3,1) = u1(1,1) + 2.5;
u1(3,2) = u1(1,2) - 4;
u1(4,1) = u1(1,1) + 8.5;
u1(4,2) = u1(1,2) - 6;
u1(5,1) = u1(1,1) + 4.5;
u1(5,2) = u1(1,2) - 8;
u1(6,1) = BS(1,1) - 2;
u1(6,2) = BS(1,2) - 20;
u1(7,1) = BS(1,1) + 2;
u1(7,2) = BS(1,2) - 18;
u1(8,1) = BS(1,1) + 2;
u1(8,2) = BS(1,2) + 22;
u1(9,1) = BS(1,1) + 6;
u1(9,2) = BS(1,2) + 14;
u1(10,1) = BS(1,1) + 8;
u1(10,2) = BS(1,2) + 22;
u1(11,1) = BS(1,1) + 24;
u1(11,2) = BS(1,2) + 5;
u1(12,1) = BS(1,1) + 26;
u1(12,2) = BS(1,2);
scatter(u1(:,1),u1(:,2),'k');
%second cell
u2 = zeros(10,2);
u2(1,1) = BS(2,1) - 20;
u2(1,2) = BS(2,2);
u2(2,1) = u2(1,1) + 8.5;
u2(2,2) = u2(1,2) + 14;
u2(3,1) = u2(1,1) + 4.5;
u2(3,2) = u2(1,2) + 8;
u2(4,1) = u2(1,1) + 9.5;
u2(4,2) = u2(1,2) - 16;
u2(5,1) = u2(1,1) + 4.5;
u2(5,2) = u2(1,2) - 8;
u2(6,1) = BS(2,1) - 2;
u2(6,2) = BS(2,2) - 20;
u2(7,1) = BS(2,1) + 6;
u2(7,2) = BS(2,2) - 14;
u2(8,1) = BS(2,1) + 2;
u2(8,2) = BS(2,2) + 22;
u2(9,1) = BS(2,1) + 6;
u2(9,2) = BS(2,2) + 8;
u2(10,1) = BS(2,1) + 16;
u2(10,2) = BS(2,2) + 19;
u2(11,1) = BS(2,1) + 24;
u2(11,2) = BS(2,2) + 5;
u2(12,1) = BS(2,1) + 18;
u2(12,2) = BS(2,2) - 8;
scatter(u2(:,1),u2(:,2),'k');
%third cell
u3(1,1) = BS(3,1) - 20;
u3(1,2) = BS(3,2);
u3(2,1) = u3(1,1) + 2.5;
u3(2,2) = u3(1,2) + 4;
u3(3,1) = u3(1,1) + 2.5;
u3(3,2) = u3(1,2) - 4;
u3(4,1) = u3(1,1) + 12.5;
u3(4,2) = u3(1,2) - 6;
u3(5,1) = u3(1,1) + 4.5;
u3(5,2) = u3(1,2) - 8;
u3(6,1) = u3(1,1) + 10.5;
u3(6,2) = u3(1,1) - 20;
u3(7,1) = u3(1,1) + 8.5;
u3(7,2) = u3(3,2) - 8;
u3(8,1) = BS(3,1) + 2;
u3(8,2) = BS(3,2) + 22;
u3(9,1) = BS(3,1) + 6;
u3(9,2) = BS(3,2) + 14;
u3(10,1) = BS(3,1) + 3;
u3(10,2) = BS(3,2) + 5;
u3(11,1) = BS(3,1) + 8;
u3(11,2) = BS(3,2) + 5;
u3(12,1) = BS(3,1) + 14;
u3(12,2) = BS(3,2) + 6;
scatter(u3(:,1),u3(:,2),'k');

% 创建 ylabel
xlabel('\fontname{宋体}\fontsize{15}横坐标\fontname{Times New Roman}\fontsize{15}/(m)',...
    'FontName','Times new roman');

% 创建 xlabel
ylabel('\fontname{宋体}\fontsize{15}纵坐标\fontname{Times New Roman}\fontsize{15}/(m)',...
    'FontName','Times new roman');

box(axes1,'on');
grid(axes1,'on');
% 设置其余坐标区属性
set(axes1,'FontName','Times new roman','FontSize',15,'YMinorTick','on');
legend('\fontname{Times new Roman}\fontsize{12}{BS1}','\fontname{Times new Roman}\fontsize{12}{BS2}','\fontname{Times new Roman}\fontsize{12}{BS3}','\fontname{Times new Roman}\fontsize{12}{User}','location','best');
axis([0 120 0 120]);
axis equal;

