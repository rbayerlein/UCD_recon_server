#!/bin/bash

path_to_output_folder=/media/rbayerlein/SSD_09_Reimund/20220218/2MLineSource_TOF_2MLineSource_TOF_171012/PET/RawData/lm_files/

path_to_lm_data=/media/rbayerlein/SSD_09_Reimund/20220218/2MLineSource_TOF_2MLineSource_TOF_171012/PET/RawData/1.2.156.112605.159303471608576.220219011013.9.6408.94037/1.2.156.112605.159303471608576.220219011352.9.13376.15564

dynamic_flag=0

framing_info=1,1799

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
