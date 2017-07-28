#include <MyOpenCal.h>
#include <OpenCAL-OMP/cal2DBuffer.h>
#include <OpenCAL-OMP/calOmpDef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

void calInit3Di(struct CALModel3D* ca3D, struct CALSubstate3Di* Q, int i, int j, int k, CALint value) {
    calSetBuffer3DElement(Q->current, ca3D->rows, ca3D->columns, i, j, k, value);
    calSetBuffer3DElement(Q->next, ca3D->rows, ca3D->columns, i, j, k, value);
}

CALint calGet3Di(struct CALModel3D* ca3D, struct CALSubstate3Di* Q, int i, int j, int k) {
    CALint ret;
    if (ca3D->is_safe == CAL_UNSAFE_ACTIVE)
        CAL_SET_CELL_LOCK_3D(i, j, k, ca3D);

    ret = calGetBuffer3DElement(Q->current, ca3D->rows, ca3D->columns, i, j, k);

    if (ca3D->is_safe == CAL_UNSAFE_ACTIVE)
        CAL_UNSET_CELL_LOCK_3D(i, j, k, ca3D);

    return ret;
}

void calInit3Dr(struct CALModel3D* ca3D, struct CALSubstate3Dr* Q, int i, int j, int k, CALreal value) {
    calSetBuffer3DElement(Q->current, ca3D->rows, ca3D->columns, i, j, k, value);
    calSetBuffer3DElement(Q->next, ca3D->rows, ca3D->columns, i, j, k, value);
}
CALreal calGet3Dr(struct CALModel3D* ca3D, struct CALSubstate3Dr* Q, int i, int j, int k) {
    CALreal ret;
    if (ca3D->is_safe == CAL_UNSAFE_ACTIVE)
        CAL_SET_CELL_LOCK_3D(i, j, k, ca3D);

    ret = calGetBuffer3DElement(Q->current, ca3D->rows, ca3D->columns, i, j, k);

    if (ca3D->is_safe == CAL_UNSAFE_ACTIVE)
        CAL_UNSET_CELL_LOCK_3D(i, j, k, ca3D);

    return ret;
}

CALint calGetX3Di(struct CALModel3D* ca3D, struct CALSubstate3Di* Q, int i, int j, int k, int n)
{
    if (ca3D->T == CAL_SPACE_FLAT)
        return calGetBuffer3DElement(Q->current, ca3D->rows, ca3D->columns, i + ca3D->X[n].i, j + ca3D->X[n].j, k + ca3D->X[n].k);
    else
        return calGetBuffer3DElement(Q->current, ca3D->rows, ca3D->columns, calGetToroidalX(i + ca3D->X[n].i, ca3D->rows), calGetToroidalX(j + ca3D->X[n].j, ca3D->columns), calGetToroidalX(k + ca3D->X[n].k, ca3D->slices));
}

CALreal calGetX3Dr(struct CALModel3D* ca3D, struct CALSubstate3Dr* Q, int i, int j, int k, int n)
{
    if (ca3D->T == CAL_SPACE_FLAT)
        return calGetBuffer3DElement(Q->current, ca3D->rows, ca3D->columns, i + ca3D->X[n].i, j + ca3D->X[n].j, k + ca3D->X[n].k);
    else
        return calGetBuffer3DElement(Q->current, ca3D->rows, ca3D->columns, calGetToroidalX(i + ca3D->X[n].i, ca3D->rows), calGetToroidalX(j + ca3D->X[n].j, ca3D->columns), calGetToroidalX(k + ca3D->X[n].k, ca3D->slices));
}

void calCopyBuffer3Dr(CALreal* M_src, CALreal* M_dest, int rows, int columns, int slices)
{
//    int tn;
//    int ttotal;
    size_t size;

//    int start;
//    int chunk;

    size = rows * columns * slices;
#ifdef GPU
#pragma omp target teams distribute parallel for num_teams(num_Teams) num_threads(num_Threads) dist_schedule(static,64) schedule(static,1)
#endif
   for(int i=0;i<size;i++)
    {
       M_dest[i]=M_src[i];
    }


//#pragma omp parallel private (start, chunk, tn, ttotal)
//    {
//        ttotal = CAL_GET_NUM_THREADS();

//        tn = CAL_GET_THREAD_NUM();
//        chunk = size / ttotal;
//        start = tn * chunk;

//        if (tn == ttotal - 1)
//            chunk = size - start;

//        memcpy(M_dest + start, M_src + start,
//               sizeof(CALreal) * chunk);
//    }
}

void calCopyBuffer3Di(CALint* M_src, CALint* M_dest, int rows, int columns, int slices)
{
//    int tn;
//    int ttotal;
    size_t size;

//    int start;
//    int chunk;

    size = rows * columns * slices;
#ifdef GPU
#pragma omp target teams distribute parallel for num_teams(num_Teams) num_threads(num_Threads) dist_schedule(static,64) schedule(static,1)
#endif
   for(int i=0;i<size;i++)
    {
       M_dest[i]=M_src[i];
    }
//#pragma omp parallel private (start, chunk, tn, ttotal)
//    {
//        ttotal = CAL_GET_NUM_THREADS();

//        tn = CAL_GET_THREAD_NUM();
//        chunk = size / ttotal;

//        start = tn * chunk;

//        if (tn == ttotal - 1)
//            chunk = size - start;

//        memcpy(M_dest + start, M_src + start,
//               sizeof(CALint) * chunk);
//    }
}

struct CALModel3D* calCADef3D(int rows,int columns,int slices, enum CALNeighborhood3D CAL_NEIGHBORHOOD_3D,enum CALSpaceBoundaryCondition CAL_TOROIDALITY,enum CALOptimization CAL_OPTIMIZATION)
{
    int i;
    struct CALModel3D *ca3D = (struct CALModel3D *)malloc(sizeof(struct CALModel3D));
    if (!ca3D)
        return NULL;

    ca3D->rows = rows;
    ca3D->columns = columns;
    ca3D->slices = slices;

    ca3D->T = CAL_TOROIDALITY;

    ca3D->OPTIMIZATION = CAL_OPTIMIZATION;
    if (ca3D->OPTIMIZATION == CAL_OPT_ACTIVE_CELLS) {
        ca3D->A.flags = calAllocBuffer3Db(ca3D->rows, ca3D->columns, ca3D->slices);
        calSetBuffer3Db(ca3D->A.flags, ca3D->rows, ca3D->columns, ca3D->slices, CAL_FALSE);
    }
    else
        ca3D->A.flags = NULL;

#pragma omp parallel
    {
#pragma omp single
             ca3D->A.num_threads = CAL_GET_NUM_THREADS();
    }


    ca3D->A.size_next = (int*)malloc(sizeof(int) * ca3D->A.num_threads);
    for(i=0;i<ca3D->A.num_threads;i++)
        ca3D->A.size_next[i] = 0;


    ca3D->A.cells = NULL;
    ca3D->A.size_current = 0;

    ca3D->X = NULL;
    ca3D->sizeof_X = 0;

    ca3D->X_id = CAL_NEIGHBORHOOD_3D;
    switch (CAL_NEIGHBORHOOD_3D) {
        case CAL_VON_NEUMANN_NEIGHBORHOOD_3D:
            calDefineVonNeumannNeighborhood3D(ca3D);
            break;
        case CAL_MOORE_NEIGHBORHOOD_3D:
            calDefineMooreNeighborhood3D(ca3D);
            break;
    }

    ca3D->pQb_array = NULL;
    ca3D->pQi_array = NULL;
    ca3D->pQr_array = NULL;
    ca3D->sizeof_pQb_array = 0;
    ca3D->sizeof_pQi_array = 0;
    ca3D->sizeof_pQr_array = 0;

    ca3D->elementary_processes = NULL;
    ca3D->num_of_elementary_processes = 0;

    ca3D->is_safe = CAL_UNSAFE_INACTIVE;

    CAL_ALLOC_LOCKS_3D(ca3D);
    CAL_INIT_LOCKS_3D(ca3D, i);


    return ca3D;
}
CALbyte* calAllocBuffer3Db(int rows, int columns, int slices) {
    return (CALbyte*)malloc(sizeof(CALbyte)*rows*columns*slices);
}
void calSetBuffer3Db(CALbyte* M, int rows, int columns, int slices, CALbyte value)
{
   memset(M, value, sizeof(CALbyte)*rows*columns*slices);
}
void calDefineVonNeumannNeighborhood3D(struct CALModel3D* ca3D)
{
    /*
         slice -1       slice 0       slice 1

       |   |         | 1 |         |   |
    ---|---|---   ---|---|---   ---|---|---
       | 5 |       2 | 0 | 3       | 6 |
    ---|---|---   ---|---|---   ---|---|---
       |   |         | 4 |         |   |
   */

    //slice  0
    calAddNeighbor3D(ca3D,   0,   0,   0);
    calAddNeighbor3D(ca3D, - 1,   0,   0);
    calAddNeighbor3D(ca3D,   0, - 1,   0);
    calAddNeighbor3D(ca3D,   0, + 1,   0);
    calAddNeighbor3D(ca3D, + 1,   0,   0);
    //slice -1
    calAddNeighbor3D(ca3D,   0,   0, - 1);
    //slice +1
    calAddNeighbor3D(ca3D,   0,   0, + 1);
}

void calDefineMooreNeighborhood3D(struct CALModel3D* ca3D)
{
    /*
         slice -1       slice 0       slice 1

    14 |10 | 17    5 | 1 | 8    23 |19 | 26
    ---|---|---   ---|---|---   ---|---|---
    11 | 9 | 12    2 | 0 | 3    20 |18 | 21
    ---|---|---   ---|---|---   ---|---|---
    15 |13 | 16    6 | 4 | 7    24 |22 | 25
   */

    //slice  0
    calAddNeighbor3D(ca3D,   0,   0,   0);
    calAddNeighbor3D(ca3D, - 1,   0,   0);
    calAddNeighbor3D(ca3D,   0, - 1,   0);
    calAddNeighbor3D(ca3D,   0, + 1,   0);
    calAddNeighbor3D(ca3D, + 1,   0,   0);
    calAddNeighbor3D(ca3D, - 1, - 1,   0);
    calAddNeighbor3D(ca3D, + 1, - 1,   0);
    calAddNeighbor3D(ca3D, + 1, + 1,   0);
    calAddNeighbor3D(ca3D, - 1, + 1,   0);
    //slice -1
    calAddNeighbor3D(ca3D,   0,   0, - 1);
    calAddNeighbor3D(ca3D, - 1,   0, - 1);
    calAddNeighbor3D(ca3D,   0, - 1, - 1);
    calAddNeighbor3D(ca3D,   0, + 1, - 1);
    calAddNeighbor3D(ca3D, + 1,   0, - 1);
    calAddNeighbor3D(ca3D, - 1, - 1, - 1);
    calAddNeighbor3D(ca3D, + 1, - 1, - 1);
    calAddNeighbor3D(ca3D, + 1, + 1, - 1);
    calAddNeighbor3D(ca3D, - 1, + 1, - 1);
    //slice +1
    calAddNeighbor3D(ca3D,   0,   0, + 1);
    calAddNeighbor3D(ca3D, - 1,   0, + 1);
    calAddNeighbor3D(ca3D,   0, - 1, + 1);
    calAddNeighbor3D(ca3D,   0, + 1, + 1);
    calAddNeighbor3D(ca3D, + 1,   0, + 1);
    calAddNeighbor3D(ca3D, - 1, - 1, + 1);
    calAddNeighbor3D(ca3D, + 1, - 1, + 1);
    calAddNeighbor3D(ca3D, + 1, + 1, + 1);
    calAddNeighbor3D(ca3D, - 1, + 1, + 1);
}


struct CALCell3D* calAddNeighbor3D(struct CALModel3D* ca3D, int i, int j, int k) {
    struct CALCell3D* X_tmp = ca3D->X;
    struct CALCell3D* X_new;
    int n;

    X_new = (struct CALCell3D*)malloc(sizeof(struct CALCell3D)*(ca3D->sizeof_X + 1));
    if (!X_new)
        return NULL;

    for (n = 0; n < ca3D->sizeof_X; n++) {
        X_new[n].i = ca3D->X[n].i;
        X_new[n].j = ca3D->X[n].j;
        X_new[n].k = ca3D->X[n].k;
    }
    X_new[ca3D->sizeof_X].i = i;
    X_new[ca3D->sizeof_X].j = j;
    X_new[ca3D->sizeof_X].k = k;

    ca3D->X = X_new;
    free(X_tmp);

    ca3D->sizeof_X++;

    return ca3D->X;
}
struct CALSubstate3Dr* calAddSubstate3Dr(struct CALModel3D* ca3D){
    struct CALSubstate3Dr* Q;
    struct CALSubstate3Dr** pQr_array_tmp = ca3D->pQr_array;
    struct CALSubstate3Dr** pQr_array_new;
    int i;

    pQr_array_new = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*(ca3D->sizeof_pQr_array + 1));
    if (!pQr_array_new)
        return NULL;

    for (i = 0; i < ca3D->sizeof_pQr_array; i++)
        pQr_array_new[i] = ca3D->pQr_array[i];

    Q = (struct CALSubstate3Dr*)malloc(sizeof(struct CALSubstate3Dr));
    if (!Q)
        return NULL;
    if (!calAllocSubstate3Dr(ca3D, Q))
        return NULL;

    pQr_array_new[ca3D->sizeof_pQr_array] = Q;
    ca3D->sizeof_pQr_array++;

    ca3D->pQr_array = pQr_array_new;
    free(pQr_array_tmp);

    return Q;
}
struct CALSubstate3Di* calAddSubstate3Di(struct CALModel3D* ca3D){
    struct CALSubstate3Di* Q;
    struct CALSubstate3Di** pQi_array_tmp = ca3D->pQi_array;
    struct CALSubstate3Di** pQi_array_new;
    int i;

    pQi_array_new = (struct CALSubstate3Di**)malloc(sizeof(struct CALSubstate3Di*)*(ca3D->sizeof_pQi_array + 1));
    if(!pQi_array_new)
        return NULL;

    for (i = 0; i < ca3D->sizeof_pQi_array; i++)
        pQi_array_new[i] = ca3D->pQi_array[i];

    Q = (struct CALSubstate3Di*)malloc(sizeof(struct CALSubstate3Di));
    if (!Q)
        return NULL;
    if (!calAllocSubstate3Di(ca3D, Q))
        return NULL;

    pQi_array_new[ca3D->sizeof_pQi_array] = Q;
    ca3D->sizeof_pQi_array++;

    ca3D->pQi_array = pQi_array_new;
    free(pQi_array_tmp);

    return Q;
}
void calInitSubstate3Di(struct CALModel3D* ca3D, struct CALSubstate3Di* Q, CALint value) {
    if (ca3D->A.cells)
    {
        calSetActiveCellsBuffer3Di(Q->current, ca3D->rows, ca3D->columns, ca3D->slices, value, ca3D->A.cells, ca3D->A.size_current);
        if(Q->next)
            calSetActiveCellsBuffer3Di(Q->next, ca3D->rows, ca3D->columns, ca3D->slices, value, ca3D->A.cells, ca3D->A.size_current);
    }
    else
    {
        calSetBuffer3Di(Q->current, ca3D->rows, ca3D->columns, ca3D->slices, value);
        if(Q->next)
            calSetBuffer3Di(Q->next, ca3D->rows, ca3D->columns, ca3D->slices, value);
    }
}
CALbyte calAllocSubstate3Dr(struct CALModel3D* ca3D,	//!< Pointer to the cellular automaton structure.
                            struct CALSubstate3Dr* Q	//!< Pointer to a 3D real (floating point) substate.
                            )
{
    Q->current = calAllocBuffer3Dr(ca3D->rows, ca3D->columns, ca3D->slices);
    Q->next = calAllocBuffer3Dr(ca3D->rows, ca3D->columns, ca3D->slices);

    if (!Q->current || !Q->next)
        return CAL_FALSE;

    return CAL_TRUE;
}


CALbyte calAllocSubstate3Di(struct CALModel3D* ca3D,struct CALSubstate3Di* Q)
{
    Q->current = calAllocBuffer3Di(ca3D->rows, ca3D->columns, ca3D->slices);
    Q->next = calAllocBuffer3Di(ca3D->rows, ca3D->columns, ca3D->slices);

    if (!Q->current || !Q->next)
        return CAL_FALSE;

    return CAL_TRUE;
}

void calSetActiveCellsBuffer3Di(CALint* M, int rows, int columns, int slices, CALint value, struct CALCell3D* active_cells, int sizeof_active_cells) {
    int n;

    for(n=0; n<sizeof_active_cells; n++)
        M[active_cells[n].k*rows*columns + active_cells[n].i*columns + active_cells[n].j] = value;
}


void calSetBuffer3Di(CALint* M, int rows, int columns, int slices, CALint value)
{
    memset(M, value, sizeof(CALint)*rows*columns*slices);
}

CALreal* calAllocBuffer3Dr(int rows, int columns, int slices) {
    return (CALreal*)malloc(sizeof(CALreal)*rows*columns*slices);
}
CALint* calAllocBuffer3Di(int rows, int columns, int slices) {
    return (CALint*)malloc(sizeof(CALint)*rows*columns*slices);
}

void calInitSubstate3Dr(struct CALModel3D* ca3D, struct CALSubstate3Dr* Q, CALreal value) {
    if (ca3D->A.cells)
    {
        calSetActiveCellsBuffer3Dr(Q->current, ca3D->rows, ca3D->columns, ca3D->slices, value, ca3D->A.cells, ca3D->A.size_current);
        if(Q->next)
            calSetActiveCellsBuffer3Dr(Q->next, ca3D->rows, ca3D->columns, ca3D->slices, value, ca3D->A.cells, ca3D->A.size_current);
    }
    else
    {
        calSetBuffer3Dr(Q->current, ca3D->rows, ca3D->columns, ca3D->slices, value);
        if(Q->next)
            calSetBuffer3Dr(Q->next, ca3D->rows, ca3D->columns, ca3D->slices, value);
    }
}

void calSetActiveCellsBuffer3Dr(CALreal* M, int rows, int columns, int slices, CALreal value, struct CALCell3D* active_cells, int sizeof_active_cells) {
    int n;

    for(n=0; n<sizeof_active_cells; n++)
        M[active_cells[n].k*rows*columns + active_cells[n].i*columns + active_cells[n].j] = value;
}

void calSetBuffer3Dr(CALreal* M, int rows, int columns, int slices, CALreal value)
{
    int size = rows * columns * slices;
    int i;

    for (i=0; i<size; i++)
        M[i] = value;
}
void calApplyElementaryProcess3D(struct CALModel3D* ca3D,CALCallbackFunc3D elementary_process)
{
    int i, j, k, n;

    if (ca3D->A.cells) //Computationally active cells optimization.
#pragma omp parallel for private (n)
        for (n = 0; n < ca3D->A.size_current; n++)
            elementary_process(ca3D, ca3D->A.cells[n].i, ca3D->A.cells[n].j, ca3D->A.cells[n].k);
    else //Standart cicle of the transition function
#pragma omp parallel for private (k, i, j)
        for (i = 0; i < ca3D->rows; i++)
            for (j = 0; j < ca3D->columns; j++)
                for (k = 0; k < ca3D->slices; k++)
                    elementary_process(ca3D, i, j, k);
}

