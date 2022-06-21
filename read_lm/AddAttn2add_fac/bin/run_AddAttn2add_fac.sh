#!/bin/bash

EXEC=/home/rbayerlein/code/explorer-master/read_lm/AddAttn2add_fac/bin/AddAttn2add_fac

lm_file=/home/rbayerlein/ssd/XPWZISNVCB
attn_dir=/home/rbayerlein/ssd/XPWZISNVCB/lm_attn_fp_exp
user_name=rbayerlein

crys_eff=/home/rbayerlein/data/explorer/20210827/Multi-Bed_Phantom_Multi-Bed_Phantom_154523/UCD/Image/crys_eff_679_840
plane_eff=/home/rbayerlein/data/explorer/20210827/Multi-Bed_Phantom_Multi-Bed_Phantom_154523/UCD/Image/plane_eff_multi_bed_phantom_tb_679x679

for (( i = 1; i < 9; i++ )); do
	sleep 1; $EXEC $lm_file/lm_reorder_f0_prompts.$i.lm $attn_dir/lm_prompts_f${i}_attn_fp_exp.raw $user_name $crys_eff $plane_eff &
done
wait

echo "$EXEC has finished."
