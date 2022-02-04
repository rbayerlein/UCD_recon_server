#!/bin/bash
echo "concatenating files..."

temp=/home/rbayerlein/ssd/MCNAVLEZST/temp
cat ${temp}/lm_reorder_f0_prompts.add_fac_out.0 ${temp}/lm_reorder_f0_prompts.add_fac_out.1 > ${temp}/lm_reorder_f0_prompts.add_fac_out


echo "done."
