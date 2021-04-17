function scanner = buildPET(name, single_ring_only)

if nargin < 2
    single_ring_only = false;
else
    if isempty(single_ring_only)
        single_ring_only = false;
    end
end

switch lower(name)
    

     case 'explorer2000mm_unitedimaging'
        name_tag = 'Explorer2000mm United-Imaging';

        ring_diameter = 786.0; % unit in mm

        crystal_size = [2.76, 2.76, 2.76, 18.1]; % tf-tr-a-d (unit in mm)
        crystal_gap_size = [0.09, 0.09, 0.09]; % tf-tr-a (unit in mm)
        
        if single_ring_only
            crystal_array_size = [35 1];
            number_of_detector_modules_axial = 1;
            detector_modula_axial_extra_offsets = [0.0];
        else
            crystal_array_size = [35 679];
            number_of_detector_modules_axial = 1;
            detector_modula_axial_extra_offsets = [0.0 0.0];
        end

        number_of_detector_modules_transaxial = 24;
        number_of_DOI_bins = 1;
        detector_module_initial_angle_offset = 0.0;
        
        number_of_projections_per_angle = 549; 

        tof_info = [430, 25];
        

    otherwise
        error('unknown scanner! micropet2, toshiba, inveon, ucdpetmr, explorer, explorer2000mm, explorer2000mm_v3_4brscanner');
        
        

        
end






% create scanner object
scanner = PETsystem(...
    name_tag, ...
    ring_diameter, ...
    crystal_size, ...
    crystal_gap_size, ...
    crystal_array_size, ...
    number_of_detector_modules_transaxial, ...
    number_of_detector_modules_axial, ...
    number_of_DOI_bins, ...
    detector_module_initial_angle_offset, ...
    detector_modula_axial_extra_offsets, ...
    number_of_projections_per_angle, ...
    tof_info);




