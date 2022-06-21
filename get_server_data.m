function get_server_data

warning off

start_dir = '/mnt/data/rbayerlein/explorer/';  


%dir_true = [handles.server_dir_backup, '/']; 
dir_true = start_dir; 
list_date_true = dir(dir_true);
dir_fake = [start_dir, 'dir_temp/']; 
mkdir(dir_fake);

for ii = 1:length(list_date_true) 
	if list_date_true(ii).isdir > 0.5 && contains(list_date_true(ii).name, '20')

		dir_date_true = [dir_true, list_date_true(ii).name, '/']; 
		dir_date_fake = [dir_fake, list_date_true(ii).name, '/'];
		mkdir(dir_date_fake); 

		list_subject_true = dir(dir_date_true); 

		for vv = 1:length(list_subject_true)
			if list_subject_true(vv).isdir > 0.5 && length(list_subject_true(vv).name) > 3
				dir_subject_true = [dir_date_true, list_subject_true(vv).name, '/']; 
				dir_subject_fake = [dir_date_fake, list_subject_true(vv).name, '/']; 
				mkdir(dir_subject_fake); 

				dir_pet_rawdata_true = [dir_subject_true, 'PET/RawData/']; 
				dir_pet_rawdata_fake = [dir_subject_fake, 'PET/RawData/']; 
				mkdir(dir_pet_rawdata_fake); 

				dir_image_true = [dir_subject_true, 'Image/']; 
				dir_image_fake = [dir_subject_fake, 'Image/']; 
				mkdir(dir_image_fake); 
				
				dir_ucd_true = [dir_subject_true, 'UCD/'];
				dir_ucd_image_true = [dir_ucd_true, 'Image/'];
				dir_ucd_pet_rawdata_true  = [dir_ucd_true, 'PET/']; 
				if ~exist(dir_ucd_true, 'dir')
					mkdir(dir_ucd_true); 	
					mkdir(dir_ucd_image_true); 
					mkdir(dir_ucd_pet_rawdata_true); 
				end
				dir_ucd_fake = [dir_subject_fake, 'UCD/']; 
				mkdir(dir_ucd_fake); 
			
				dir_ucd_image_fake = [dir_subject_fake,'UCD/Image/']; 
				mkdir(dir_ucd_image_fake); 

				dir_ucd_pet_rawdata_fake = [dir_subject_fake, 'UCD/PET/']; 
				mkdir(dir_ucd_pet_rawdata_fake); 
				
				

				list_pet_rawdata_true = dir(dir_pet_rawdata_true); 

				for jj = 1:length(list_pet_rawdata_true) 
					if list_pet_rawdata_true(jj).isdir > 0.5 && length(list_pet_rawdata_true(jj).name) > 3
						dir_pet_acq_true = [dir_pet_rawdata_true, list_pet_rawdata_true(jj).name,'/']; 
						dir_pet_acq_fake = [dir_pet_rawdata_fake, list_pet_rawdata_true(jj).name,'/'];
						mkdir(dir_pet_acq_fake); 

						lm_list = dir(dir_pet_acq_true); 

						for qq = 1:length(lm_list)
							if contains(lm_list(qq).name, '.1.dcm')
								dcm_copy = [dir_pet_acq_true, lm_list(qq).name]; 
								copyfile(dcm_copy, dir_pet_acq_fake); 
							end
							if contains(lm_list(qq).name, '.1.raw')
								raw_create = [dir_pet_acq_fake, lm_list(qq).name];
								fid_1 = fopen(raw_create, 'w'); 
								fwrite(fid_1, 1, 'float'); 
								fclose(fid_1); 
							end
						end
					end
				end


				list_image_true = dir(dir_image_true); 
				for uu = 1:length(list_image_true)
					if list_image_true(uu).isdir > 0.5 && contains(list_image_true(uu).name, 'CT') && contains(list_image_true(uu).name, 'AC')
						dir_image_create = [dir_image_fake, list_image_true(uu).name,'/']; 
						mkdir(dir_image_create); 
					end
					%if contains(list_image_true(uu).name,'.sen_img')
						%sen_create = [dir_image_fake, list_image_true(uu).name]; 
						%fid_2 = fopen(sen_create, 'w'); 
						%fwrite(fid_2, 1, 'float'); 
						%fclose(fid_2); 
					%end
				end


				list_ucd_image_true = dir(dir_ucd_image_true); 
				for uu = 1:length(list_ucd_image_true)
					if contains(list_ucd_image_true(uu).name,'.sen_img')
						sen_create_ucd = [dir_ucd_image_fake, list_ucd_image_true(uu).name]; 
						fid_3 = fopen(sen_create_ucd, 'w'); 
						fwrite(fid_3, 1, 'float'); 
						fclose(fid_3); 
					end
				end
				

				
			end
		end
	end
end





quit; 

