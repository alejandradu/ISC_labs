%function [vterm_high] = high_re (r, pO, pF, vis)
%    [vterm_high] = (2*g^2*r^2*(pO - pF))/(9*vis);
%end

%function [vterm_low] = low_re (r, pO, pF)
%    [vterm_low] = (8*g*r*(pO - pF)/(3*cd_high*pF))^(1/2)
%end

l_bac = 10^-1;
v_bac = 10^-6;
l_hum = 1.7;
v_hum = 1.5;
n_water = 0.1002;
p_water = 998;
n_gly = 1.412;
p_gly = 1260;
p_al = 2790;
p_st = 8000;
g = 9.81;
cd_high = 0.4;

diameters = [1/16 5/64 3/32 1/8];
r = diameters * (2.54/200); %in meters

%pO is aluminun and pF is glycerol
%vterm is in m/s
%converted

vterm_low = (2*g^2*r.^2*(p_al - p_gly))/(9*n_gly)*100;

%pO is steel and pF is water. Cd is 0.4 for high Re
%vterm is in m/s
%converted 

vterm_high= (8*g*(p_st - p_water)/(3*cd_high*p_water))^0.5*(r.^0.5);
vterm_high = vterm_high*100;

radius = [0.00005:0.0001:0.0025]; %in meters

vlow_model = (2*g^2*(radius.^2)*(p_al - p_gly))/(9*n_gly)*100;
vhigh_model = (8*g*(p_st - p_water)/(3*cd_high*p_water))^0.5*(radius.^0.5);
vhigh_model = vhigh_model*100;

r = r*100;
radius = radius*100;

figure('Name', 'terminal velocity vs. radius size. Low Re')

subplot(1,2,1)
hold on
plot(radius, vlow_model, 'LineWidth', 1.5)
scatter(r, vterm_low, 100, [0.5843    0.8157    0.9882], 'filled')
hold off
xlabel('Radius (cm)')
ylabel('Terminal velocity (cm/s)')
title('Low Reynolds environment')
set(gca,'FontSize',19)
set(gcf,'color','w');

subplot(1,2,2)
hold on
plot(radius, vhigh_model, 'LineWidth', 1.5)
scatter(r, vterm_high, 100, [0.5843    0.8157    0.9882], 'filled')
hold off
xlabel('Radius (cm)')
ylabel('Terminal velocity (cm/s)')
title('High Reynolds environment')
set(gca,'FontSize',19)
set(gcf,'color','w');