%A class to perform fast listmode- or sinogram-based PET image reconstruction
%
%Author: Jian Zhou
%
%Revised: Xuezhu Zhang
%
%Date: 2013
%



classdef PETsystem
    
    properties(GetAccess = 'public', SetAccess = 'private')
        name_tag; % scanner name
        system_parms; % structure that contains various system parameters
    end
    
    
    
    
    
    
    
    methods(Access = 'public')
        
        
        
        
        function obj = PETsystem(...
                tag, ...
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
                tof_info)
            
            obj.name_tag = tag;
            obj.system_parms.ring_diameter = ring_diameter;
            obj.system_parms.crystal_size = crystal_size;
            obj.system_parms.crystal_gap_size = crystal_gap_size;
            obj.system_parms.crystal_array_size = crystal_array_size;
            obj.system_parms.number_of_detector_modules_transaxial = ...
                number_of_detector_modules_transaxial;
            obj.system_parms.number_of_detector_modules_axial = ...
                number_of_detector_modules_axial;
            obj.system_parms.number_of_DOI_bins = number_of_DOI_bins;
            obj.system_parms.detector_module_initial_angle_offset = ...
                detector_module_initial_angle_offset;
            obj.system_parms.detector_module_axial_extra_offsets = ...
                detector_modula_axial_extra_offsets;
            obj.system_parms.number_of_projections_per_angle = ...
                number_of_projections_per_angle;
            obj.system_parms.tof_info = tof_info;
            obj.system_parms.projector = 'Siddon';
            obj.system_parms.depth_ratio = 0.5;
        end
        
        
        
        
        
        function obj = setProjector(obj, projector)
            %set projector type: `Siddon' or `Bresenham' or 'Linterp'
            %function obj = setProjector(projector)
            obj.system_parms.projector = projector;
        end
        
        
        
        
        
        
        function obj = setDepthRatio(obj, ratio)
            %set a ratio to determine average depth interaction point inside crystal
            %0.5 is default, meaning at the crystal center
            %function obj = setDepthRatio(ratio)
            if ratio < 0 || ratio > 1.0
                error('ratio must be in [0 1]!');
            end
            obj.system_parms.depth_ratio = ratio;
        end
        
        
        
        
        
        function obj = setNumberOfProjectionsPerAngle(obj, num_of_radial_bins)
            %set number of radial bins
            %function obj = setNumberOfProjectionsPerAngle(num_of_radial_bins)
            obj.system_parms.number_of_projections_per_angle = num_of_radial_bins;
        end
        
        
        
        
        
        function dmo = getDetectorModuleAxialOffsets(obj)
            %calculate detector module axial offsets
            %function dmo = getDetectorModuleAxialOffsets
            crystal_axial_pitch = obj.system_parms.crystal_size(3) + ...
                obj.system_parms.crystal_gap_size(3);
            detector_module_axial_size = crystal_axial_pitch * ...
                obj.system_parms.crystal_array_size(2);
            
            nb = obj.system_parms.number_of_detector_modules_axial;
            dmo = (- nb * 0.5 + 0.5 + (0 : (nb-1))) * ...
                detector_module_axial_size + ...
                obj.system_parms.detector_module_axial_extra_offsets;
        end
        
        
        
        
        
        
        function nr = getNumberOfCrystalRings(obj)
            %get number of crystal rings
            %function nr = getNumberOfCrystalRings
            nr = obj.system_parms.crystal_array_size(2) * ...
                obj.system_parms.number_of_detector_modules_axial;
        end
        
        
        
        
        
        
        function na = getDefaultNumberOfAngles(obj)
            %get number of projection angles in default mode
            %function na = getDefaultNumberOfAngles
            na = obj.system_parms.crystal_array_size(1) * ...
                obj.system_parms.number_of_detector_modules_transaxial / 2;
        end
        
        
        
        
        
        
        function ro = getCrystalRingOffsets(obj)   % define axial gap (2014-0311)
            %get crystal ring axial offsets
            %function ro = getCrystalRingOffsets
            crystal_axial_pitch = obj.system_parms.crystal_size(3) + ...
                obj.system_parms.crystal_gap_size(3);
            dm_offsets = getDetectorModuleAxialOffsets(obj);
            nb = obj.system_parms.crystal_array_size(2);
            xt_centers = (- nb * 0.5 + 0.5 + (0 : (nb-1))) * crystal_axial_pitch;
            ro = [];
            for n = 1 : obj.system_parms.number_of_detector_modules_axial
                ro = [ro , xt_centers + dm_offsets(n)];
            end
        end
        
        
        
        
        
        
        function tc = getCrystalTransaxialLocations(obj)
            %calculate crystal bin transaxial coordinates (bin means DOI bin)
            %function tc = getCrystalTransaxialLocations
            tfs = (obj.system_parms.crystal_size(1) + obj.system_parms.crystal_gap_size(1));
            trs = (obj.system_parms.crystal_size(2) + obj.system_parms.crystal_gap_size(2));
            nxtal_trans = obj.system_parms.crystal_array_size(1);
            df = (-nxtal_trans*0.5 + 0.5 + (0 : nxtal_trans-1)) * tfs;
            dr = (-nxtal_trans*0.5 + 0.5 + (0 : nxtal_trans-1)) * trs;
            
            num_of_doi_bins = obj.system_parms.number_of_DOI_bins;
            xtal_loc = zeros(2, num_of_doi_bins, nxtal_trans);
            xtal_size_depth = obj.system_parms.crystal_size(4);
            for n=1: nxtal_trans
                l = df(n)-dr(n);
                dl = l / num_of_doi_bins;
                cl = l * 0.5 - ((0:num_of_doi_bins-1) + 0.5) * dl;
                h = xtal_size_depth;
                dh = h / num_of_doi_bins;
                ch = h * 0.5 - ((0:num_of_doi_bins-1) + 0.5) * dh;
                xtal_loc(:, :, n) = [cl(:) + dr(n), ch(:)]';
            end
            
            R = obj.system_parms.ring_diameter * 0.5;
            nblock_trans = obj.system_parms.number_of_detector_modules_transaxial;
            tc = zeros(2, num_of_doi_bins, nxtal_trans, nblock_trans);
            for n=1:nblock_trans
                % always start at 3 o'clock
                a0 = (n-1) * 2*pi / nblock_trans  + ...
                    obj.system_parms.detector_module_initial_angle_offset * pi / 180;
                xs = xtal_loc(1,:,:);
                ys = xtal_loc(2,:,:);
                
                t = pi * 0.5 + a0;
                x = xs * cos(t) - ys * sin(t) + (R + xtal_size_depth*obj.system_parms.depth_ratio) * cos(a0);
                y = xs * sin(t) + ys * cos(t) + (R + xtal_size_depth*obj.system_parms.depth_ratio) * sin(a0);
                
                tc(1,:,:,n) = x;
                tc(2,:,:,n) = y;
            end
            tc = reshape(tc, 2, num_of_doi_bins, nxtal_trans*nblock_trans);
            tc = squeeze(tc);
        end
        
        
        
        
        
        
        
        function plot(obj, options)
            % plot scanner geometry with options: 2d,trans,simple; 2d,trans,simple,+lor; 3d,simple
            %function plot(options)
            switch lower(options)
                case '2d,geom'
                    line_color = [0, 0, 0];
                    c_width_front = obj.system_parms.crystal_size(1) + obj.system_parms.crystal_gap_size(1);
                    c_width_back = obj.system_parms.crystal_size(2) + obj.system_parms.crystal_gap_size(2);
                    b_width_front = obj.system_parms.crystal_array_size(1) * c_width_front;
                    b_width_back = obj.system_parms.crystal_array_size(1) * ...
                        (obj.system_parms.crystal_size(2) + obj.system_parms.crystal_gap_size(2));
                    b_depth = obj.system_parms.crystal_size(end);
                    for n = 1 : (obj.system_parms.crystal_array_size(1)+1)
                        pfront(n,:) = [-b_depth/2, b_width_front/2 - (n-1)*c_width_front];
                        pback(n,:) = [b_depth/2, b_width_back/2 - (n-1)*c_width_back];
                    end
                    
                    nblock_trans = obj.system_parms.number_of_detector_modules_transaxial;
                    a0 = obj.system_parms.detector_module_initial_angle_offset;
                    cc = obj.system_parms.ring_diameter * 0.5 + b_depth / 2;
                    plot(0, 0, '+'); hold on;
                    for m = 1 : nblock_trans
                        
                        a = 360 / nblock_trans * pi / 180 * (m - 1) + a0;
                        r = [cos(a), -sin(a); sin(a), cos(a)];
                        t = [cc * cos(a); cc * sin(a)];
                        
                        pf = r * pfront' + repmat(t, 1, size(pfront,1));
                        pb = r * pback' + repmat(t, 1, size(pback, 1));
                        for n = [1:(obj.system_parms.crystal_array_size(1)+1)]
                            line([pf(1,n), pb(1,n)], [pf(2,n), pb(2,n)], 'color', line_color);
                        end
                        line([pf(1,1), pf(1,end)], [pf(2,1), pf(2,end)], 'color', line_color);
                        line([pb(1,1), pb(1,end)], [pb(2,1), pb(2,end)], 'color', line_color);
                        
                    end
                    maxc = ceil((cc+b_depth)/100)*100;
                    line([-maxc,maxc], [0, 0], 'linestyle', '--', 'color', 'k');
                    line([0, 0], [-maxc,maxc], 'linestyle', '--', 'color', 'k');
                    hold off;
                    axis equal;
                    set(gca, 'drawMode', 'fast', ...
                        'xtickmode', 'manual', ...
                        'ytickmode', 'manual');
                    xlim([-maxc,maxc]);
                    ylim([-maxc,maxc]);
                    set(gca, 'xtick', [-maxc,0,maxc]);
                    set(gca, 'ytick', [-maxc,0,maxc]);
                    xlabel('x (mm)', 'fontsize', 12, 'fontweight', 'bold');
                    ylabel('y (mm)', 'fontsize', 12, 'fontweight', 'bold');
                    set(gca, 'fontsize', 12);
                    title(sprintf('%s - 2D Geometry', obj.name_tag), 'fontweight', 'bold');
                    %            		grid on;
                    %            		set(gca, 'gridlinestyle', '--');
                    
                case '2d,trans,simple'
                    tc = getCrystalTransaxialLocations(obj);
                    if size(tc, 3) > 1
                        for n = 1 : size(tc, 2)
                            plot(squeeze(tc(1,n,:)), squeeze(tc(2,n,:)), '.');
                            if n==1
                                hold on;
                            end
                        end
                        
                        R = obj.system_parms.ring_diameter * 0.5 - 25;
                        nblock_trans = obj.system_parms.number_of_detector_modules_transaxial;
                        for n = 1 : nblock_trans
                            a = 360/nblock_trans * (n-1) + ...
                                obj.system_parms.detector_module_initial_angle_offset;
                            text(R*cos(a*pi/180), R*sin(a*pi/180), num2str(n), ...
                                'horizontalalignment', 'center');
                        end
                        axis square;
                        
                        hold off;
                    else
                        plot(tc(1,:), tc(2,:), '.');
                        R = obj.system_parms.ring_diameter * 0.5 - 25;
                        nblock_trans = obj.system_parms.number_of_detector_modules_transaxial;
                        for n = 1 : nblock_trans
                            a = 360/nblock_trans * (n-1) + ...
                                obj.system_parms.detector_module_initial_angle_offset;
                            text(R*cos(a*pi/180), R*sin(a*pi/180), num2str(n), ...
                                'horizontalalignment', 'center');
                        end
                        title(sprintf('%s - 2D Geometry', obj.name_tag), 'fontweight', 'bold');
                        axis square;
                    end
                    
                case '2d,trans,simple,+lor'
                    tc = getCrystalTransaxialLocations(obj);
                    if size(tc, 3) > 1
                        for n = 1 : size(tc, 2)
                            plot(squeeze(tc(1,n,:)), squeeze(tc(2,n,:)), '.');
                            if n==1
                                hold on;
                            end
                        end
                        
                        axis square;
                        hold off;
                    else
                        plot(tc(1,:), tc(2,:), '.');
                        axis square;
                    end
                    xp = getDefaultSinogramCrystalPairs(obj);
                    np = obj.system_parms.number_of_projections_per_angle;
                    xp = reshape(xp, 2, np, size(xp,2)/np);
                    hold on;
                    for i = 1 : floor(size(xp,3)) : size(xp, 3)
                        for n = 1 : np
                            xtal_pair = xp(:, n, i);
                            p0 = tc(:, xtal_pair(1));
                            p1 = tc(:, xtal_pair(2));
                            line([p0(1), p1(1)], [p0(2), p1(2)]);
                        end
                    end
                    hold off;
                    title('Only show LORs at the 1st angle');


                    


                case '2d,trans,simple_gap,+lor'


                    crystal_array = obj.system_parms.crystal_array_size(1)
                    num_blocks = obj.system_parms.number_of_detector_modules_transaxial
                    
                    num_crystals = crystal_array * num_blocks

                    % crystal_no = 1:num_crystals;

                    crystal_no_filled = [];

                    crystal_no_filled_offset = [];

                    zeros_gaps = zeros(1, crystal_array);


                    for no_block_filled = 1:num_blocks

                        crystal_array_no_block_filled = 1 + (no_block_filled-1)*crystal_array : crystal_array + (no_block_filled-1)*crystal_array;

                        if mod(no_block_filled, 2) ~= 0
                            crystal_no_filled = [crystal_no_filled, crystal_array_no_block_filled];
                        else
                            crystal_no_filled = [crystal_no_filled, zeros_gaps];
                        end

                        if mod(no_block_filled, 2) == 0
                            crystal_no_filled_offset = [crystal_no_filled_offset, crystal_array_no_block_filled];
                        else
                            crystal_no_filled_offset = [crystal_no_filled_offset, zeros_gaps];
                        end

                    end

                    whos crystal_no_filled crystal_no_filled_offset



                    % xp = int16(obj.getDefaultSinogramCrystalPairs());
                    xp = getDefaultSinogramCrystalPairs(obj);

                    xp_1 = size(xp, 1);
                    xp_2 = size(xp, 2);

                    xp_oneblockgap = zeros(xp_1, xp_2, 'int16');
                    xp_oneblockgap_offset = zeros(xp_1, xp_2, 'int16');


                    for nx = 1:xp_1
                        for ny = 1:xp_2

                            xp_oneblockgap(nx, ny) = crystal_no_filled(xp(nx, ny));
                            xp_oneblockgap_offset(nx, ny) = crystal_no_filled_offset(xp(nx, ny));

                        end
                    end


                    tc = getCrystalTransaxialLocations(obj);
                    if size(tc, 3) > 1
                        for n = 1 : size(tc, 2)
                            plot(squeeze(tc(1,n,:)), squeeze(tc(2,n,:)), '.');
                            if n==1
                                hold on;
                            end
                        end
                        
                        axis square;
                        hold off;
                    else
                        plot(tc(1,:), tc(2,:), '.');
                        axis square;
                    end


                    % xp = getDefaultSinogramCrystalPairs(obj);
                    % xp = xp_oneblockgap;
                    xp = xp_oneblockgap_offset;

                    np = obj.system_parms.number_of_projections_per_angle;
                    xp = reshape(xp, 2, np, size(xp,2)/np);
                    % whos xp
                    % max(xp(:))
                    % min(xp(:))

                    % np3 = size(xp, 3);

                    % xp = xp(xp>0);

                    % whos xp

                    % np2 = size(xp, 2);

                    % xp = reshape(xp, 1, np2, length(xp,1)/np2);


                    hold on;
                    for i = 1 : floor(size(xp,3)) : size(xp, 3)
                        for n = 1 : np  %  np2  % np
                            xtal_pair = xp(:, n, i);

                            xtal_pair1 = xtal_pair(1);
                            xtal_pair2 = xtal_pair(2);

                            % whos xtal_pair

                            if xtal_pair(1)>0 & xtal_pair(2)>0

                                % xtal_pair = xtal_pair(xtal_pair(1)>0, xtal_pair(2)>0);
                                % whos xtal_pair

                                p0 = tc(:, xtal_pair(1));
                                p1 = tc(:, xtal_pair(2));
                                line([p0(1), p1(1)], [p0(2), p1(2)]);

                            end

                        end
                    end
                    hold off;
                    title('Only show LORs at the 1st angle');
                    

                    
                case '3d,simple'
                    tc = getCrystalTransaxialLocations(obj);
                    ro = getCrystalRingOffsets(obj);
                    for z = 1 : length(ro)
                        if size(tc, 3) > 1
                            for n = 1 : size(tc, 2)
                                plot3(squeeze(tc(1,n,:)), ...
                                    squeeze(tc(2,n,:)), ...
                                    ro(z) * ones(size(tc,3),1), '.');
                                if n==1 && z==1
                                    hold on;
                                end
                            end
                            axis square;
                        else
                            plot3(tc(1,:), tc(2,:), ro(z) * ones(size(tc,2),1), '.');
                            axis square;
                            if z==1
                                hold on;
                            end
                        end
                    end
                    hold off;
                otherwise
                    error('unknown options! valid options: "2d,geom", "2d,trans,simple", "2d,trans,simple,+lor", "3d,simple"');
            end
        end
        
        
        
        
        
        
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    methods
        
        
        
        
        
        function xtal_pairs = getDefaultSinogramCrystalPairs(obj)
            % create crystal pairs according my own numbering scheme
            %function xtal_pairs = getDefaultCrystalPairs
            %
            
            % always equal to total number of crystals per ring divided by
            % 2
            nxtal_trans = obj.system_parms.crystal_array_size(1);
            xtal_num_trans_total = obj.system_parms.number_of_detector_modules_transaxial * ...
                obj.system_parms.crystal_array_size(1);
            
            num_of_angles = xtal_num_trans_total / 2;
            
            if obj.system_parms.number_of_projections_per_angle > (fix(xtal_num_trans_total/2)*2 - 2)
                error('too many projections!');
            end
            
            % case when detector module is not located exactly at 3-clock
            if obj.system_parms.detector_module_initial_angle_offset ~= 0
                disp('The 1st detector module is rotated by an angle!');
                disp('NOTE: For crystal pairing, this angle offset is never used directly!');
                disp('but always assume it is equal to 180 / (number_of_detector_modules_transaxial)');
                id0 = 0;
                id1 = num_of_angles;
                nr = fix(num_of_angles/2) * 2 - 2;
                h0 = zeros(nr,1);
                h1 = zeros(nr,1);
                
                for i=0:nr-1
                    
                    id=id0 + fix(nr/2) - i - 1;
                    if id < 0
                        id = id + xtal_num_trans_total;
                    end
                    
                    h0(i+1)=id;
                    id = id1-fix(nr/2) + i;
                    h1(i+1)=id;
                end
            else
                
                id0 = fix(nxtal_trans / 2);
                odd = mod(nxtal_trans, 2) ~= 0;
                
                if odd
                    id1 = fix(xtal_num_trans_total / 2) + fix(nxtal_trans / 2);
                    nr = fix(num_of_angles / 2) * 2 - 1;
                else
                    id1 = fix(xtal_num_trans_total / 2) + fix(nxtal_trans / 2) - 1;
                    nr = fix(num_of_angles / 2) * 2;
                end
                
                h0 = zeros(nr,1);
                h1 = zeros(nr,1);
                
                for i=0:nr-1
                    if odd
                        id = id0 + fix(nr/2) - i;
                    else
                        id = id0 + fix(nr/2) -1 - i;
                    end
                    
                    if id < 0
                        id = id + xtal_num_trans_total;
                    end
                    
                    h0(i+1) = id;
                    if odd
                        id = id1 - fix(nr/2) + i;
                    else
                        id = id1 - fix(nr/2) + 1 + i;
                    end
                    h1(i+1) = id;
                end
            end
            
            % pairing
            c = 1; k = 1;
            while c < length(h0)
                xtal_id1 = h0(c);
                xtal_id2 = h1(c);
                xp_first_angle(:,k) = [xtal_id1; xtal_id2];
                k = k + 1;
                if (c+1) <= length(h1)
                    xtal_id1 = h0(c);
                    xtal_id2 = h1(c+1);
                    xp_first_angle(:,k) = [xtal_id1; xtal_id2];
                    k = k + 1;
                end
                c = c + 1;
            end
            
            %
            nn = size(xp_first_angle,2);
            num_of_projs_per_angle = obj.system_parms.number_of_projections_per_angle;
            xtal_pairs = zeros(2, num_of_angles * num_of_projs_per_angle);
            k = 1;
            for i=1:num_of_angles
                for j=1:num_of_projs_per_angle
                    pp = xp_first_angle(:, fix(nn / 2) - fix(num_of_projs_per_angle / 2) + j);
                    p0 = mod(pp(1) + i-1, xtal_num_trans_total);
                    p1 = mod(pp(2) + i-1, xtal_num_trans_total);
                    xtal_pairs(:,k) = [p0; p1];
                    k = k + 1;
                end
            end
            xtal_pairs = xtal_pairs + 1;
        end
        
        
        
        
        
        
        
        
        
        
        
        %         function block_pairs_sinogram = getDefaultSinogramBlockPairs(obj, nblockwidth, nblockwidth_half)
        %
        %
        %             block_num_trans_total = obj.system_parms.number_of_detector_modules_transaxial;
        %
        %             num_of_angles = block_num_trans_total / 2;
        %
        %
        %             for i = 1:num_of_angles
        %
        %
        %
        %
        %             end
        %
        %         end
        
        
        
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    methods(Access = 'public')
        
        
        
        function rp = getDefaultSinogramPlaneOrder(obj)
            %get default plane arrangement in terms of ring pair (different from Micxx-gram)
            %function rp = getDefaultSinogramPlaneOrder
            %
            nring = obj.getNumberOfCrystalRings();
            rp = zeros(2, nring*nring);
            offset = 0;
            for n = 1 : nring
                if n==1
                    rp(:,1:nring) = [1:nring; 1:nring];
                    offset = offset + nring;
                else
                    r_odd = [1:(nring-n+1); n:(nring)];
                    r_even = nring-r_odd+1;
                    nr = (nring-n+1)*2;
                    rp(:, offset + (1:2:nr)) = r_odd;
                    rp(:, offset + (2:2:nr)) = r_even;
                    offset = offset + nr;
                end
            end
        end
    end
    
    
    
    
    
    
    
    methods(Access = 'public')
        
        
        
        function prjs = doListModeForwardProjectionNonTOF(obj, image, image_size, voxel_size, lmdata)
            %do listmode-based forward projection for non-TOF case
            %function prjs = doListModeForwardProjectionNonTOF(image, image_size, voxel_size, lmdata)
            %
            image = reshape(image, image_size);
            %
            if (~isa(lmdata, 'int16')) && (~isa(lmdata, 'uint16'))
                error('invalid data type: lmdata must be int16 or uint16!');
            end
            %
            tc = obj.getCrystalTransaxialLocations();
            to = obj.getCrystalRingOffsets();
            switch obj.system_parms.projector
                case 'Siddon'
                    tic
                    prjs = fproj_mt(image, image_size, voxel_size, tc, to, lmdata);
                    toc
                case 'Bresenham'
                    tic
                    prjs = fproj_mt_bresenham(image, image_size, voxel_size, tc, to, lmdata);
                    toc
                case 'Linterp'
                    tic
                    prjs = fproj_mt_linterp(image, image_size, voxel_size, tc, to, lmdata);
                    toc
                otherwise
                    error('invalid projector!');
            end
        end
        


       
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
end
