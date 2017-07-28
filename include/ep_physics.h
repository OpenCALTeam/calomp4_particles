#ifndef EP_PHYSICS_H
#define EP_PHYSICS_H

#include <model.h>
#pragma omp declare target
void applyForce(CALreal* F, CALreal* p0, CALreal* v0, CALreal m, CALreal t, CALreal* pf, CALreal* vf);
void resetF(CALint* ID_current,CALreal* Q_next,int cell_x, int cell_y, int cell_z);
#pragma omp end declare target
#endif
