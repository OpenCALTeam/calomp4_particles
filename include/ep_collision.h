#ifndef EP_COLLISION_H
#define EP_COLLISION_H

#include <model.h>
#pragma omp declare target
void inner_collision(CALint* ID_current,CALreal* Q_current,CALreal* Q_next, int cell_x, int cell_y, int cell_z);
void outer_collision(CALint* ID_current,CALreal* Q_current,CALreal* Q_next,CALint* Xi,CALint* Xj,CALint* Xk, int cell_x, int cell_y, int cell_z);
#pragma omp end declare target
#endif
