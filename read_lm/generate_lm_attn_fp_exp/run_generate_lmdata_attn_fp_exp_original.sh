#!/bin/bash

FIL=/home/raid12/zhxie/EMIC-phantom/2021-04-20/Davis_recon/lmrecon_exploer/app/lm_fp_exp

if [ -f $FIL ]; then
    	echo " Compiled Recon program  exists in $FIL." 
else  
    	echo "Compiled Recon program doesn't exists in $FIL, please specify the correct directory." 
fi
 
#for (( m=1; m<=8; m++ )) do
for (( m=8; m>=1; m-- )) do
{
	FIL=./lm_attn_fp_exp/lm_prompts_f$[$m]_attn_fp_exp.raw
	
	if [ -f $FIL ]; then
    	echo "File $FIL exists."    	
	else
	/home/raid12/zhxie/EMIC-phantom/2021-04-20/Davis_recon/lmrecon_exploer/app/lm_fp_exp   ./cfg_attn_fp_exp/lmrecon_attn_fp_exp_f$m.cfg

	fi
}
done


