#!/bin/bash

P1=\'/mnt/data/rbayerlein/data/motion/20191212/COOPER_JEANEE_1841349_091702//Image/CTAC_201.sen_img\'
P2=\'/mnt/data/rbayerlein/data/motion/20191212/COOPER_JEANEE_1841349_091702//Image//CTAC_201_mumap_kVp-140_size-256x256x646_vox-2.7344x2.7344x3.img\'
P3=\'/mnt/data/rbayerlein/data/motion/20191212/COOPER_JEANEE_1841349_091702//Image/crys_eff_679x840\'
P4=\'/mnt/data/rbayerlein/data/motion/20191212/COOPER_JEANEE_1841349_091702//Image/plane_eff_679x679\'

matlab -nodesktop -nojvm -r "cd /mnt/data/rbayerlein/code/explorer-master/reconstruction/lmrecon_senimg/; make_senimg2(${P1},${P2},${P3},${P4}); quit"
