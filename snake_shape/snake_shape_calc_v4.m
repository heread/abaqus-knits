clear
close all

%v1: WORKING!!!
%v2: r + L INDEPENDENT!!
%v3: can do 3 rows (prob not n rows)
%v4: gonna try to do n rows

%todo: get this working for arbitrary L (ugh)

r = 1/(2*pi);
char_L = 2*pi*r; %THIS IS T
% whole_char_L = 15; %4 * pi * r is the "natural length" of one bend
whole_char_L = 100;

total_L_set = whole_char_L;
num_rows = 3;

a = @(p) 3*pi/2 - p; %p \in [-pi/2, pi/2] --> [2pi, pi]
g = @(t,p) t*(a(p) - 2*pi) + 2*pi; 
%p = -pi/2: [0,1] --> [2pi, 2pi]
%p = pi/2: [0,1] --> [2pi, pi]

x = @(t,p) r*(p + g(t,p) - sin(g(t,p)));
y = @(t,p) r*(1 - cos(g(t,p)));

T = 2*pi*r; %period for one cycloid
circ_angles = [-pi/2, pi/2];
bd_mid = T + r*circ_angles; %[left bd, right bd]
bd_end = 2*T;

t_slice = 0:0.01:1; %corr. to time in g(t) 0 --> 1
t_sliceCoarse = 0:0.1:1;

n_p = @(t,p) r*(t*cos(t*(p + pi/2)) - t + 1);
m_p = @(t,p) r*t*sin(t*(p + pi/2));


input_all = linspace(0,total_L_set,30);
x_all = 0:0.01:total_L_set;
% x_all = 0:0.01:12; %todo: change this back
all_dist = zeros(1,length(t_slice));

figure(1)
hold on

% BLUE CURVES
for i = 1:length(input_all)
    ret_cur = calc_curves(t_slice,input_all(i),r,total_L_set,num_rows);
    plot(ret_cur(1,:),ret_cur(2,:), 'b--')
end

%PINK CURVES
for i = 1:length(t_sliceCoarse)
    ret_cur = calc_curves(t_sliceCoarse(i),x_all,r,total_L_set,num_rows);
    %plot(ret_cur(1,:),ret_cur(2,:),'m','Linewidth',1.5)
    plot(ret_cur(1,:),ret_cur(2,:),'Linewidth',1.5)
    all_dist(i) = calc_dist(ret_cur(1,:),ret_cur(2,:));
end

%TESTINNG CURVES
% for i = 1:length(t_slice)
%     ret_cur = calc_curves(t_slice(i),x_all,r,total_L_set,num_rows);
%     all_dist(i) = calc_dist(ret_cur(1,:),ret_cur(2,:));
% end

% 
% ylim([-0.01,2])
% xlim([-0.01,2])
xlim([-1,6])
axis equal

% figure(2)
% plot(t_slice,all_dist)
% xlabel('t')
% ylabel('curve length')