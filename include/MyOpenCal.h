#ifndef MYOPENCAL_H
#define MYOPENCAL_H

#include <OpenCAL-OMP/cal3D.h>
#include <OpenCAL-OMP/cal3DIO.h>
#include <omp.h>

struct CALModel3D* calCADef3D(int rows,int columns,int slices,enum CALNeighborhood3D CAL_NEIGHBORHOOD_3D,enum CALSpaceBoundaryCondition CAL_TOROIDALITY,enum CALOptimization CAL_OPTIMIZATION);
#pragma omp declare target
#define calSetBuffer3DElement(M, rows, columns, i, j, k, value) ( (M)[( ((k)*(rows)*(columns)) + ((i)*(columns)) + (j) )] = (value) )
#define calGetBuffer3DElement(M, rows, columns, i, j, k) ( M[( ((k)*(rows)*(columns)) + ((i)*(columns)) + (j) )] )

CALreal calGet3Dr(struct CALModel3D* ca3D, struct CALSubstate3Dr* Q, int i, int j, int k);
CALint calGet3Di(struct CALModel3D* ca3D, struct CALSubstate3Di* Q, int i, int j, int k);
void calInit3Di(struct CALModel3D* ca3D, struct CALSubstate3Di* Q, int i, int j, int k, CALint value);
void calInit3Dr(struct CALModel3D* ca3D, struct CALSubstate3Dr* Q, int i, int j, int k, CALreal value);
CALint calGetX3Di(struct CALModel3D* ca3D, struct CALSubstate3Di* Q, int i, int j, int k, int n);
CALreal calGetX3Dr(struct CALModel3D* ca3D, struct CALSubstate3Dr* Q, int i, int j, int k, int n);
void calCopyBuffer3Dr(CALreal* M_src, CALreal* M_dest, int rows, int columns, int slices);
void calCopyBuffer3Di(CALint* M_src, CALint* M_dest, int rows, int columns, int slices);
CALbyte* calAllocBuffer3Db(int rows, int columns, int slices);
void calSetBuffer3Db(CALbyte* M, int rows, int columns, int slices, CALbyte value);
void calDefineVonNeumannNeighborhood3D(struct CALModel3D* ca3D);
void calDefineMooreNeighborhood3D(struct CALModel3D* ca3D);
struct CALCell3D* calAddNeighbor3D(struct CALModel3D* ca3D, int i, int j, int k);
struct CALSubstate3Dr* calAddSubstate3Dr(struct CALModel3D* ca3D);
struct CALSubstate3Di* calAddSubstate3Di(struct CALModel3D* ca3D);
void calInitSubstate3Di(struct CALModel3D* ca3D, struct CALSubstate3Di* Q, CALint value);
CALbyte calAllocSubstate3Dr(struct CALModel3D* ca3D,struct CALSubstate3Dr* Q);
struct CALSubstate3Di* calAddSubstate3Di(struct CALModel3D* ca3D);
CALbyte calAllocSubstate3Di(struct CALModel3D* ca3D,struct CALSubstate3Di* Q);
void calSetActiveCellsBuffer3Di(CALint* M, int rows, int columns, int slices, CALint value, struct CALCell3D* active_cells, int sizeof_active_cells);
void calSetBuffer3Di(CALint* M, int rows, int columns, int slices, CALint value);
CALreal* calAllocBuffer3Dr(int rows, int columns, int slices);
CALint* calAllocBuffer3Di(int rows, int columns, int slices);
void calInitSubstate3Dr(struct CALModel3D* ca3D, struct CALSubstate3Dr* Q, CALreal value);
void calSetActiveCellsBuffer3Dr(CALreal* M, int rows, int columns, int slices, CALreal value, struct CALCell3D* active_cells, int sizeof_active_cells);
void calSetBuffer3Dr(CALreal* M, int rows, int columns, int slices, CALreal value);
void calApplyElementaryProcess3D(struct CALModel3D* ca3D,CALCallbackFunc3D elementary_process);


#pragma omp end declare target




#endif
