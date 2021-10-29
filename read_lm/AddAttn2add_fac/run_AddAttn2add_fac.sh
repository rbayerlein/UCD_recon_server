#!/bin/bash

EXEC=AddAttn2add_fac

lm_file=/home/rbayerlein/ssd/ADJTOHCQKP
attn_dir=/home/rbayerlein/ssd/ADJTOHCQKP/
user_name=rbayerlein
nc_file=/home/rbayerlein/data/explorer/20200124/cylinder_phantom3_E-20200124-172216-11_172342/PET/RawData/1.2.156.112605.18587648329783.200125012343.9.6520.125134/1.2.156.112605.18587648329783.200125012731.9.13476.18310.1.nc
crys_eff=/home/rbayerlein/data/explorer/20200124/cylinder_phantom3_E-20200124-172216-11_172342/UCD/Image/crys_eff_679_840
plane_eff=/home/rbayerlein/data/explorer/20200124/cylinder_phantom3_E-20200124-172216-11_172342/UCD/Image/plane_eff_679x679

for (( i = 1; i < 9; i++ )); do
	sleep 1; ./$EXEC $lm_file/lm_reorder_f0_prompts.$i.lm $attn_dir/lm_prompts_f${i}_attn_fp_exp.raw $user_name $crys_eff $plane_eff &
done
wait

echo "$EXEC has finished."
