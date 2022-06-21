function cod_xyz = get_ROI_COD(handles)

p1 = handles.r1_pos_all(handles.num_roi,:); 
p2 = handles.r2_pos_all(handles.num_roi,:);
    
x_span = [floor(p1(2)), floor(p1(2)+p1(4))];
y_span = [floor(p2(2))-256, floor(p2(2)+p2(4)) - 256];
z_span = [floor(p1(1)), floor(p1(1)+p1(3))];

roi_img = zeros(size(handles.mask_img)); 
mask_img(y_span(1):y_span(2), x_span(1):x_span(2), z_span(1):z_span(2)) = 1; 

roi_img = roi_img .* handles.mask_img; 
roi_img(roi_img < 0) = 0; 

roi_list = roi_img(roi_img > 0.5); 
roi_list = roi_list(:); 


rois = unique(roi_list); 

wts = ones(size(rois)); 
for k = 1:length(rois)
	r_temp = roi_list == rois(k); 
	wts(k) = length(r_temp); 
end

wts = wts ./ sum(wts(:)); 



for k = 1:length(rois)
	cod_temp = handles.cod_xyz(:, rois(k):handles.num_patch:end); 

	cod_xyz = cod_xyz + cod_temp .* wts(k); 
end



