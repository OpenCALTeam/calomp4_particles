#ifndef EP_MOVILI_CAZZU_H
#define EP_MOVILI_CAZZU_H

#include <model.h>
#pragma omp declare target
void moviliCazzu(CALint* ID_current,CALint* ID_next,CALreal* Q_current,CALreal* Q_next,CALint* Xi,CALint* Xj,CALint* Xk, int cell_x, int cell_y, int cell_z);
void pezziala(CALint* ID_next,CALreal* Q_next,int slot,int cell_x, int cell_y, int cell_z);
void sucala(CALint* ID_current,CALint* ID_next,CALreal* Q_current,CALreal* Q_next,CALint* Xi,CALint* Xj,CALint* Xk,int destination_slot, int source_slot, int cell_x, int cell_y, int cell_z, int n);
#pragma omp end declare target
#endif
