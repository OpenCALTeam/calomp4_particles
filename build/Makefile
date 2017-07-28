CC=/home/pnapoli/clang/bin/clang++
CC_RELEASE=/home/pnapoli/clang_release/bin/clang++ 
WARNING= -w
OPENMP=-fopenmp=libomp
CUDA=-fopenmp-targets=nvptx64-nvidia-cuda
SOURCE=../source_cli/cli.c ../source/*.c 
HEADER= -I ../include/
OPENCALHEADER=-I /usr/local/opencal-1.0/include/
OPENCALSOURCE=/home/pnapoli/Desktop/opencal-1.0/OpenCAL-OMP/source/*.c 
OUTPUT_D= -o cal_dem_d
OUTPUT_R= -o cal_dem_r
OUTPUT_CPU_P_R= -o cal_dem_c_p_r
OUTPUT_CPU_R= -o cal_dem_c_r

#$(OPENMP) $(CUDA) $(WARNING)

cpu_parallel_release:
	$(CC_RELEASE) $(SOURCE) $(HEADER) $(OPENCALHEADER) $(OPENMP) $(OUTPUT_CPU_P_R)
cpu_release:
	$(CC_RELEASE) $(SOURCE) $(HEADER) $(OPENCALHEADER) $(OUTPUT_CPU_R)
debug: 
	$(CC) $(SOURCE) $(HEADER) $(OPENCALHEADER) $(OPENMP) $(CUDA) $(OUTPUT_D) 

release:
	
	$(CC_RELEASE) $(SOURCE) $(HEADER) $(OPENCALHEADER) $(OPENMP) $(CUDA) $(OUTPUT_R) 

list:
	@echo "List of commands of this MakeFile\n"
	@echo "CPU commands:"
	@echo "\tcpu_parallel_release"
	@echo "\tcpu_release"
	@echo "\nGPU commands:"
	@echo "\tdebug"
	@echo "\trelease\n"
	@echo "\nUtility:"
	@echo "\tclean"


clean:
	rm cal_dem_*

