#!/bin/bash

path_to_output_folder=/mnt/ssd/rbayerlein/explorer/recon_data_temp/phantom_6beds_78rings_from188_34rings_overlap_220_4it_scatter_corr_20210916/

path_to_lm_data=/mnt/data/rbayerlein/explorer/20210827/Multi-Bed_Phantom_Multi-Bed_Phantom_154523/PET/RawData/1.2.156.112605.159303471608576.210827224523.9.6628.91186/1.2.156.112605.159303471608576.210827225240.9.12756.14170

dynamic_flag=0

framing_info=1,1320

write_lm_flag=1

write_3d_sino_flag=0

write_blocksino_flag=0

write_tof_histo_flag=0

for i in {1..8..1}; do
	echo $path_to_output_folder > Reconstruction_Parameters_$i
	echo $path_to_lm_data.$i.raw >> Reconstruction_Parameters_$i
	echo $dynamic_flag >> Reconstruction_Parameters_$i
	echo $framing_info >> Reconstruction_Parameters_$i
	echo $write_lm_flag >> Reconstruction_Parameters_$i
	echo $write_3d_sino_flag >> Reconstruction_Parameters_$i
 	echo $write_blocksino_flag >> Reconstruction_Parameters_$i
	echo $write_tof_histo_flag >> Reconstruction_Parameters_$i
done
