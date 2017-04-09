%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% introduction 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Purpose: MATLAB script to calculate the settlement of soil around each tunnel section
assumptions:
	> centre-to-centre distance (c-t-c d) between running tunnels are 30m (TFL, 2013)
	> centre-to-centre distance (c-t-c d) between tunnels at platforms are 30m (TFL, 2013)
	> settlement from indiviudal tunnels are added together with superposition
	> for each section, tunnel have the same radius and depth
	> all soils are fine-grained (check)
	> worst case scenario between O'Reilly and New Method & practical method was chosen
		worst case conditions are determined by: slope
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a csv file with following information for each section to design
% X-coord, Y-coord, depth, radius, (c-t-c d) b/w tunnels, K (practical purpose), Vs/Vt ratio
tunnel_sections = dlmread('tunnel sections.csv',',',0,0);

%% face loss
Vl = [0.75, 1.5];			% unit: % -> for assignment use 0.75 and 1.5

%% accuracy
increment_size = 0.01;
round_decimal = 2;

%% settlement points to find
settlement_points = 0:0.0025:0.04;
settlement_slope_lim = 1/500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iterate each section and calculate the settlement curve and contour
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,~] = size(tunnel_sections);
settlement_contour_plot_Vl = {};
for section_id = 1:m
	tunnels_sel = tunnel_sections(section_id, :);
    settlement_curve_plot_Vl = {};
	for face_loss = 1:2
		%% empty matrix for collectin of data
		settlement_curve_Vl = [];
		settlement_contour_Vl = [];
		other_calc_values = []; 	% V0, Vl, Vt, Vs, i_rn, i_pr, Smax
		settlement_slope_Vl = [];

		%% V
		V0 = pi*tunnels_sel(4)^2;
		Vt = V0*Vl(face_loss)/100;
		Vs = Vt*tunnels_sel(7);

		%% i
		% O'Reilly and New Method
		i_rn = 0.43*tunnels_sel(3) + 1.1;
		% practical purpose
		i_pr = tunnels_sel(3)*tunnels_sel(6);
		% worst case - determined by smaller i value
		i_cr = min(i_rn,i_pr);

		%% Smax
		Smax = Vs/(sqrt(2*pi)*i_cr);

		%% tabulate calculated values
		other_calc_values = [V0, Vl(face_loss), Vt, Vs, i_rn, i_pr, Smax];

		%% settlement curve
        y_cr = [-3*i_cr,(ceil(-(3*i_cr)/increment_size)*increment_size):increment_size:(floor((3*i_cr)/increment_size)*increment_size),3*i_cr];
		Sy_L_cr = [];
		Sy_R_cr = [];
		Sy_total_cr = [];
		
		for aa = 1:length(y_cr)
			Sy_L_cr(aa,:) = [y_cr(aa)-(tunnels_sel(5)/2), Smax*exp(-(y_cr(aa)^2)/(2*i_cr^2))];
			Sy_R_cr(aa,:) = [y_cr(aa)+(tunnels_sel(5)/2), Smax*exp(-(y_cr(aa)^2)/(2*i_cr^2))];
		end

		yy_cr = round([ceil(min(Sy_L_cr(:,1))/increment_size)*increment_size:increment_size:floor(max(Sy_R_cr(:,1))/increment_size)*increment_size],round_decimal,'decimals');
        Sy_total_cr(1,:) = [min(Sy_L_cr(:,1)), 0];
        for bb = 1:length(yy_cr)
			% Left tunnel
			if isempty(find(round(Sy_L_cr(:,1),round_decimal,'decimals')==yy_cr(bb),1)) == 0
				Sy_LL_cr = Sy_L_cr(find(round(Sy_L_cr(:,1),round_decimal,'decimals')==yy_cr(bb),1),2);
			else
				Sy_LL_cr = 0;
			end

			% Right tunnel
			if isempty(find(round(Sy_R_cr(:,1),round_decimal,'decimals')==yy_cr(bb),1)) == 0
				Sy_RR_cr = Sy_R_cr(find(round(Sy_R_cr(:,1),round_decimal,'decimals')==yy_cr(bb),1),2);
			else
				Sy_RR_cr = 0;
			end

			% total
			Sy_total_cr(bb+1,:) = [yy_cr(bb), round(Sy_LL_cr+Sy_RR_cr,4,'decimals')];
        end
        Sy_total_cr(length(Sy_total_cr)+1,:) = [max(Sy_R_cr(:,1)), 0];
        
		settlement_curve_Vl = Sy_total_cr;

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% X and Y coordinates for contour lines
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% train lines
        if section_id == m
        	m_tunnel = (tunnels_sel(2) - tunnel_sections(section_id-1, 2))/(tunnels_sel(1) - tunnel_sections(section_id-1, 1));
        else
        	m_tunnel = (tunnels_sel(2) - tunnel_sections(section_id+1, 2))/(tunnels_sel(1) - tunnel_sections(section_id+1, 1));
        end

        SY = Sy_total_cr(:,1);
		SSy = Sy_total_cr(:,2);
        k = 1;
		for sett_point = 1:length(settlement_points)
			Sy = find(SSy == settlement_points(sett_point));
			if isempty(Sy) == 0
				for z = 1:length(Sy)
					x_cont = tunnels_sel(1) - SY(Sy(z))/sqrt(1+(1/m_tunnel.^2));
					y_cont = (-x_cont/m_tunnel)+(tunnels_sel(2)+(tunnels_sel(1)/m_tunnel));
                    settlement_contour(k,:) = [SY(Sy(z)), round(SSy(Sy(z)),4,'decimals'), x_cont, y_cont];
                    k = k+1;
                end
			end
        end

		settlement_contour_Vl_sort = settlement_contour;
	    [~,order] = sort(settlement_contour_Vl_sort(:,1));        
	    settlement_contour_Vl = settlement_contour_Vl_sort(order,:);

	    settlement_contour_plot_Vl{section_id,face_loss} = settlement_contour_Vl;

	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%% differential settlement
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		Y = settlement_contour_Vl(:,1);
		ss = settlement_contour_Vl(:,2);

		L = [0];
		delta_L = [0];
        L_delta = [0];
        check = [0];
		for diff_sett_points = 2:length(Y)
			L(diff_sett_points) = abs(Y(diff_sett_points)-Y(diff_sett_points-1));
			delta_L(diff_sett_points) = abs(ss(diff_sett_points)-ss(diff_sett_points-1))/L(diff_sett_points);
            L_delta(diff_sett_points) = 1/delta_L(diff_sett_points);
            if delta_L(diff_sett_points) < settlement_slope_lim
                check(diff_sett_points) = 1;
            else
                check(diff_sett_points) = 0;
            end
        end

		settlement_slope_Vl = [Y, ss, L', delta_L', L_delta', check'];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% export data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculated values
        dlmwrite(strcat('calculated_values_',int2str(section_id),'_',int2str(face_loss),'.csv'),other_calc_values,'delimiter',',');
        dlmwrite(strcat('settlement_curve_',int2str(section_id),'_',int2str(face_loss),'.csv'),settlement_curve_Vl,'delimiter',',');
        dlmwrite(strcat('settlement_contour_',int2str(section_id),'_',int2str(face_loss),'.csv'),settlement_contour_Vl,'delimiter',',');
        dlmwrite(strcat('settlement_slope_',int2str(section_id),'_',int2str(face_loss),'.csv'),settlement_slope_Vl,'delimiter',',');
        
        settlement_curve_plot_Vl{face_loss} = [settlement_curve_Vl];
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% plot the graph
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	face_1 = settlement_curve_plot_Vl{1};
% 	face_2 = settlement_curve_plot_Vl{2};
% 
% 	figure(section_id)
% 	ax = axes; 
% 	ax.YDir = 'reverse';
% 	hold on
% 	plot(face_1(:,1),face_1(:,2),'r-')
% 	plot(face_2(:,1),face_2(:,2),'b--')
% 	legend('face loss = 0.75%','face loss = 1.5%','Location','southwest')
% 	xlabel('y[m]')
% 	ylabel('settlement[m]')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% filter the contour coordinates into settlements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for face_losses = 1:2
	for sett_points = 1:length(settlement_points)
		contour_XY = [];
        kk = 1;
		for section_ids = 1:m
			concerned_section = settlement_contour_plot_Vl{section_ids, face_losses};
			Sy_cont = find(concerned_section(:,2) == settlement_points(sett_points));
			if isempty(Sy_cont) == 0
				for zz = 1:length(Sy_cont)
                    contour_XY(kk,:) = [section_ids, settlement_points(sett_points), concerned_section(Sy_cont(zz),3), concerned_section(Sy_cont(zz),4)];
                    kk = kk+1;
                end
			end
		end
		dlmwrite(strcat('settlement_contours_XY_',int2str(sett_points),'_',int2str(face_losses),'.csv'),contour_XY,'delimiter',',');
	end
end

disp('done')