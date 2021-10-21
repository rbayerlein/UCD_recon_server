#!/bin/bash

EXEC=/home/rbayerlein/code/explorer-master/read_lm/lmrecon_exploer/app/lm_fp_exp

if [ -f $EXEC ]; then
    	echo " Compiled Recon program  exists in $EXEC." 
else  
    	echo "Compiled Recon program doesn't exists in $EXEC, please specify the correct directory." 
    	exit 1
fi
 
#for (( m=1; m<=8; m++ )) do
for (( m=8; m>=1; m-- )) do
{
	FIL=./lm_attn_fp_exp/lm_prompts_f$[$m]_attn_fp_exp.raw
	
	if [ -f $FIL ]; then
    	echo "Attenuation lm file $FIL exists already."    	
	else
		echo "Attenuation lm file $FIL does not exist yet. Will create now..."
		$EXEC ./cfg_attn_fp_exp/lmrecon_attn_fp_exp_f$m.cfg
	fi
}
done
# check if all files exist (the empty ones are not created by default)
for (( m=8; m>=1; m-- )) do
{
	FIL=./lm_attn_fp_exp/lm_prompts_f$[$m]_attn_fp_exp.raw
	if [ -f $FIL ]; then
		echo "Attenuation lm file $FIL successfully created."
	else
		echo "Attenuation lm file $FIL was not created (probably because the original lm file was empty). Will create dummy file."
		touch $FIL
	fi
}
echo "run_generate_lmdata_attn_fp_exp.sh done."
done