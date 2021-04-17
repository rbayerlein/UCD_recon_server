CC=g++
INC=../inc
SRC_PATH=../src

SRC=$(SRC_PATH)/petsys_cfg.cpp $(SRC_PATH)/petsys_cscanner.cpp $(SRC_PATH)/petsys_lmprj.cpp $(SRC_PATH)/petsys_ipsf.cpp $(SRC_PATH)/petsys_spmtx.cpp $(SRC_PATH)/petsys_reg.cpp $(SRC_PATH)/petsys_log.cpp $(SRC_PATH)/petsys_prj.cpp $(SRC_PATH)/petsys_timer.cpp $(SRC_PATH)/petsys_recon.cpp


SRC1=$(SRC_PATH)/petsys_cfg.cpp $(SRC_PATH)/petsys_cscanner.cpp $(SRC_PATH)/petsys_lmprj_any_noaddmul.cpp $(SRC_PATH)/petsys_ipsf.cpp $(SRC_PATH)/petsys_spmtx.cpp $(SRC_PATH)/petsys_reg.cpp $(SRC_PATH)/petsys_log.cpp $(SRC_PATH)/petsys_prj.cpp $(SRC_PATH)/petsys_timer.cpp $(SRC_PATH)/petsys_recon.cpp


lmrecon_nontof:
	$(CC) -o lmrecon_nontof lmrecon.cpp $(SRC) -I$(INC) -O3 -DDEBUG -fopenmp -DUSE_OMP -lm   # -DUSE_TOF 

lmrecon_tof:
	$(CC) -o lmrecon_tof lmrecon.cpp $(SRC) -I$(INC) -O3 -DUSE_TOF -DDEBUG -fopenmp -DUSE_OMP -lm


