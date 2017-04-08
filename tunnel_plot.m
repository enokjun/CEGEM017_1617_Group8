%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tunnel
tunnel = 'train line 1';
r = 4.3; 		% radius of excavated circular area  unit:m

% Full slip
Sn_fs = 10.66;		% unit: kpa
St_fs = 0;          % unit: kpa

% No slip
Sn_ns = -11.84;     % unit: kpa
St_ns = 45;         % unit: kpa

% Nu effect
Nu = 746.84;		% unit: kN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% angle
% angle theta = 0 starts from the top-most point
% the angle increase in clock-wise direction 
theta_max = 90;
theta_increment = 10;

% size
linethickenss = 2;
tunnel_linethickness = 5;
title_font_size_scale_factor = 1;
legend_font_size = 10;

% r limit
rmin_M = -40;
rmax_M = 40;

rmin_N = 600;
rmax_N = 900;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% angles matrix 
theta = 0:theta_max;

% full slip condition
N_fs = ((-r/3)*(Sn_fs+2*St_fs)*cosd(2*theta))+Nu;	% force along the tunnel lining (unit: kN/m)
M_fs = ((r*r)/6)*(2*Sn_fs+St_fs)*cosd(2*theta);	% momemt along the tunnel lining (unit: kNm/m)

% no slip condition
N_ns = ((-r/3)*(Sn_ns+2*St_ns)*cosd(2*theta))+Nu;	% force along the tunnel lining (unit: kN/m)
M_ns = ((r*r)/6)*(2*Sn_ns+St_ns)*cosd(2*theta);	% momemt along the tunnel lining (unit: kNm/m) 

% Nu
Nu_mat = repmat(Nu,length(theta),1);

% Tunnel Lining Location
tunnel_lining = zeros(length(theta),1);

%% plot N - full slip, no slip and Nu
figure(1)
pax = polaraxes;
pax.Title.String = strcat(tunnel,{' : '},'Compression Force [kN/m of tunnel]');
pax.TitleFontSizeMultiplier = title_font_size_scale_factor;
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
pax.ThetaLim = [0 theta_max];
pax.ThetaTick = [0:theta_increment:theta_max];
pax.RLim = [rmin_N rmax_N];
hold on
polarplot(deg2rad(theta),N_fs,'b-','LineWidth',linethickenss)
polarplot(deg2rad(theta),N_ns,'r--','LineWidth',linethickenss)
polarplot(deg2rad(theta),Nu_mat,'k-','LineWidth',linethickenss)
legend({'full slip','no slip','Nu effect'},'Fontsize',legend_font_size)
% rlabel('Force [kN/m]')
% thetalabel('Angle [degree]')

%% plot M - full slip and no slip
% r limit
%{
rmin = round(min(min(M_fs),min(M_ns)),1,'significant');
rmax = round(max(max(M_fs),max(M_ns)),1,'significant');
Nofdigit_rmin = 1+floor(log10(abs(rmin)));
Nofdigit_rmax = 1+floor(log10(abs(rmax)));
rmin = floor(rmin/(10^Nofdigit_rmin))*(10^Nofdigit_rmin);
rmax = ceil(rmax/(10^Nofdigit_rmax))*(10^Nofdigit_rmax);
%}
figure(2)
pax2 = polaraxes;
pax2.Title.String = strcat(tunnel,{' : '},'Bending Moment [kNm/m of tunnel]');
pax2.TitleFontSizeMultiplier = title_font_size_scale_factor;
pax2.ThetaZeroLocation = 'top';
pax2.ThetaDir = 'clockwise';
pax2.ThetaLim = [0 theta_max];
pax2.ThetaTick = [0:theta_increment:theta_max];
pax2.RLim = [rmin_M rmax_M];
hold on
polarplot(deg2rad(theta),M_fs,'b-','LineWidth',linethickenss)
polarplot(deg2rad(theta),M_ns,'r--','LineWidth',linethickenss)
polarplot(deg2rad(theta),tunnel_lining,'k-','LineWidth',tunnel_linethickness)
legend({'full slip','no slip','tunnel lining'},'Fontsize',legend_font_size)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% max and min 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('maximum N and M for full slip:')
disp(max(N_fs))
disp(max(M_fs))
disp('minimum N and M for full slip:')
disp(min(N_fs))
disp(min(M_fs))

disp('maximum N and M for no slip:')
disp(max(N_ns))
disp(max(M_ns))
disp('minimum N and M for no slip:')
disp(min(N_ns))
disp(min(M_ns))