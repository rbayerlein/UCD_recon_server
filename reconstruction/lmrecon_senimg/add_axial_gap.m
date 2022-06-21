function [crys_eff_wgap, plane_eff_wgap] = add_axial_gap(crys_eff, plane_eff)

num_crys_ax_wogap = 672; 
num_crys_ax_wgap = 679; 
num_crys_trans = 840;
num_crys_ax_unit = 84;
num_unit = 8;  

size_1 = size(crys_eff)
size_2 = size(plane_eff)

if (size_1(1) + size_1(2)) ~= (num_crys_ax_wogap + num_crys_trans)
	disp('reshaping crystal efficiency...');
	crys_eff = reshape(crys_eff, num_crys_ax_wogap, num_crys_trans); 
end

if (size_2(1) + size_2(2)) ~= (num_crys_ax_wogap + num_crys_ax_wogap)
	disp('reshaping plane efficiency...')
	plane_eff = reshape(plane_eff, num_crys_ax_wogap, num_crys_ax_wogap); 
end

crys_eff_wgap = zeros(num_crys_ax_wgap, num_crys_trans); 
counter = 1; 
for k = 0:(num_unit-1)
	crys_eff_wgap((counter+k):(counter+k+num_crys_ax_unit-1),:) = crys_eff(counter:(counter+num_crys_ax_unit-1),:); 
	counter = counter + num_crys_ax_unit; 
end


plane_eff_wgap = zeros(num_crys_ax_wgap, num_crys_ax_wgap); 
for k = 0:(num_crys_ax_wogap-1)
	for kk = 0:(num_crys_ax_wogap-1)
		ax_new1 = floor(k/num_crys_ax_unit) + k; 
		ax_new2 = floor(kk/num_crys_ax_unit) + kk; 
		plane_eff_wgap(ax_new1+1, ax_new2+1) = plane_eff(k+1, kk+1); 
	end
end

figure
imagesc(crys_eff_wgap)

figure
imagesc(plane_eff_wgap)
