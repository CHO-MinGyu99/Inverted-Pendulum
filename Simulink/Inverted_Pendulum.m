clc;
clear;
close all;
syms s f;
%% Parameter change
% M = 1 * 0.9;
% m = 0.3 * 0.9;
% g = 981 * 0.9;
% l = 30 * 0.9;
% b = 0.5 * 0.9;
% r = 1.27 * 0.9;
% km = 4.9 * 0.9;
% kb = 0.0507 * 0.9;
% R = 0.3 * 0.9;
% 
% I = (m*l^2)/3;
% alpha = (M+m)*(I+m*l^2)-(m*l)^2;
% Fr = b+(2*pi/r)^2*(km*kb/R);
% Fv = (2*pi/r)*(km/r);
% 
% F = [0      1              0           0
%      0 -(I+m*l^2)*Fr*m*l/alpha  (m^2*g*l^2)/alpha   0
%      0      0              0           1
%      0 -(m*l*Fr)/alpha       m*g*l*(M+m)/alpha  0];
% G = [     0
%      (I+m*l^2)*Fv/alpha
%           0
%         m*l*Fr/alpha];
% H = [1 0 0 0
%      0 0 1 0];
% J = [0
%      0];
%% state-space
F = [0     1      0    0
     0 -19.2817 52.693 0
     0     0      0    1
     0 -0.2709  14.525 0];
G = [   0
     75.0201
        0
      1.0541];
H = [1 0 0 0
     0 0 1 0];
J = [0
     0];

H_cart = H(1,:);
H_angle = H(2,:);

states = {'x' 'x_dot' 'phi' 'phi_dot'};
inputs = {'u'};
outputs = {'x'; 'phi'};

sys = ss(F,G,H,J,'statename',states,'inputname',inputs,'outputname',outputs); % 상태공간모델 생성

poles = eig(F);             % pole
transFc = tf(sys);          % system

zero_Group = transFc.Numerator;
pole_Group = transFc.Denominator;

x_trans_num = conv(zero_Group{1,1},zero_Group{1,1});
x_trans_den = conv(pole_Group{1,1},pole_Group{1,1});
sys_x = tf(x_trans_num,x_trans_den);

angle_trans_num = conv(zero_Group{2,1},zero_Group{2,1});
angle_trans_den = conv(pole_Group{2,1},pole_Group{2,1});
sys_angle = tf(angle_trans_num,angle_trans_den);
%% Controllability
C = [G F*G (F^2)*G (F^3)*G];
determinent_C = det(C);
co      = ctrb(sys); % 가제어성행렬 생성    
rank_co = rank(co);     % Rank를 통한 가제어성 판별
%% Control canonical form
invC = inv(C);
tn = [0 0 0 1] * invC;
invTc = [tn*F*F*F; tn*F*F; tn*F; tn];
Tc = inv(invTc);
Ac = invTc * F * Tc;
Bc = invTc * G;
Cc = H * Tc;
Dc = J;
%% Stabilizing Controller
% Not satisfied
% p = [-2.5556-2.5356i -2.5556+2.5356i -25.56 -38.334];
% s1 = p(1);
% s2 = p(2);
% s3 = p(3);
% s4 = p(4);
% a_c = (f-p(1))*(f-p(2))*(f-p(3))*(f-p(4));
% a_c = sym2poly(a_c);
% a_cF = F^4 + a_c(2)*F^3 + a_c(3)*F^2 + a_c(4)*F + a_c(5)*eye(4,4);
% K = [0 0 0 1] * invC * a_cF;

% satisfied
pole = [-2.4-1.8i -2.4+1.8i -24 -36];
s1 = pole(1);
s2 = pole(2);
s3 = pole(3);
s4 = pole(4);
K = acker(F,G,[s1 s2 s3 s4]);
% bode(ss(F-G*K,G,H,0));
%% Reference Controller
ref = 5; % 5번에서는 20사용 / 6번에서는 5사용
Hn = H_cart;  % only cart Position
Jn = [0];
N = inv([F G; Hn Jn]) * [0;0;0;0;1];
Nx = N(1:4);
Nu = N(5);
Nbar = Nu + K * Nx;
[num,den] = ss2tf(F-G*K,G,Hn,Jn);
sys_cl = tf(num,den);
% pzmap(sys_cl);
eig(F-G*K);
stepinfo(ref*sys_cl);
%% SRL(LQR) Controller
% x에 대한 SRL

% rlocus(sys_x);
% [rho_x LQR_x_pole]=rlocfind(sys_x)
rho_x = 0.8188;
LQR_pole = [-19.8737+3.3166i -19.8737-3.3166i -3.7028+0.0117i -3.7028-0.0117i];
K_LQR = acker(F,G,LQR_pole);
[num_LQR,den_LQR] = ss2tf(F-G*K_LQR,G,Hn,0);
sys_LQR = tf(num_LQR,den_LQR);
Nbar_LQR = Nu + K_LQR * Nx;
stepinfo(sys_LQR);

% Q, R을 선정해 lqr함수를 이용한 LQR

Q = [30  0  0  0
     0  0  0  0
     0  0  300000  0
     0  0  0  0];
R = 1;
%K_LQR = lqr(F,G,Q,R);
%% Observability
%O = [H; H*F; H*F^2; H*F^3];
O = obsv(F,H);
rank_O = rank(O);
%% Estimator
p_e = pole * 10;
L = place(F',H',p_e)';
%% Estimator SRL
p_e_LQR = LQR_pole * 10;
L_LQR = place(F',H',p_e_LQR)';
%% Integral Controller
FI = [0 Hn;zeros(4,1) F];
GI = [0; G];
K_temp_I = acker(FI,GI,[LQR_pole(1) LQR_pole(2) LQR_pole(3) LQR_pole(4) LQR_pole(4)]);
K_Integral = [K_temp_I(2) K_temp_I(3) K_temp_I(4) K_temp_I(5)];
KI_integral = K_temp_I(1);
%% Best Combination
K_temp_I_SRLX = acker(FI,GI,[pole(1) pole(2) pole(3) pole(4) pole(4)]);
K_Integral_SRLX = [K_temp_I_SRLX(2) K_temp_I_SRLX(3) K_temp_I_SRLX(4) K_temp_I_SRLX(5)];
KI_integral_SRLX = K_temp_I_SRLX(1);