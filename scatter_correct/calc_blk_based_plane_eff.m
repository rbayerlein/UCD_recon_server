function fname_plane_eff_file_blk = calc_blk_based_plane_eff(plane_eff_file)
% function to calcualte the plane efficiency on block level
% used for normalization of block based randoms sinograms

%% fixed parameters
num_crys_ax_w_gap = 679;
num_crys_ax_wo_gap = 672;
num_cry_per_unit = 84;
num_ax_crys_per_block = 6;
num_ax_blks = 112;  %8 units * 14 blocks per unit

%% open file
plaeff_crys = fread(fopen(plane_eff_file, 'rb'), inf, 'float');
plaeff_crys = reshape(plaeff_crys, num_crys_ax_w_gap, num_crys_ax_w_gap);

%% main 
plaeff_blks = zeros(num_ax_blks, num_ax_blks);

for axA = 1 : num_crys_ax_w_gap
    for axB = 1 : num_crys_ax_w_gap
        if rem(axA,num_cry_per_unit+1) == 0 || rem(axB,num_cry_per_unit+1) == 0
            continue;
        end
        unitA = ceil(axA/num_cry_per_unit);
        unitB = ceil(axB/num_cry_per_unit);
        crysA = axA-unitA+1;
        crysB = axB-unitB+1;
        blkA=ceil(crysA/num_ax_crys_per_block);
        blkB=ceil(crysB/num_ax_crys_per_block);
        plaeff_blks(blkA, blkB) = plaeff_blks(blkA, blkB) + plaeff_crys(axA,axB)/(power(num_ax_crys_per_block,2));
    end
end
% figure;
% imshow(plaeff_blks, []);

%% saving
dim = strfind(plane_eff_file, '679x679');
fname_plane_eff_file_blk = [plane_eff_file(1:dim-1), 'blk_112x112'];
fid_out = fopen(fname_plane_eff_file_blk, 'wb');
fwrite(fid_out, plaeff_blks, 'float');
fclose(fid_out);

end % function