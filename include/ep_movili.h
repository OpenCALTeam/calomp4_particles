#ifndef EP_MOVILI_H
#define EP_MOVILI_H

#include <model.h>
#pragma omp declare target
void movili(CALint* ID_current,CALreal* Q_current,CALreal* Q_next, int cell_x, int cell_y, int cell_z);
#pragma omp end declare target
#endif
