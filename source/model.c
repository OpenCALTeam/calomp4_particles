#include <boundary.h>
#include <ep_collision.h>
#include <ep_movili.h>
#include <ep_movili_cazzu.h>
#include <ep_physics.h>
#include <ep_utils.h>
#include <init.h>
#include <model.h>
#include <utils_io.h>
#include <sim_stop.h>
#include <stdlib.h>

struct CALModel3D* u_modellu;
struct Substates Q;

CALreal* Q_current = NULL;
CALreal* Q_next = NULL;
CALint* ID_current = NULL;
CALint* ID_next = NULL;

CALint* Xi = NULL;
CALint* Xj = NULL;
CALint* Xk = NULL;
CALint step;
CALint initial_nummber_of_particles;
CALreal elapsed_time;
#define num_Threads 64
#define num_Teams 128
//#define VERBOSE
#define GPU
//------------------------------------------------------------------

void mapperToSubstates3D(struct CALModel3D *model, CALreal * realSubstate_current_OUT, CALint* intSubstate_current_OUT) {

  int ssNum_r = REAL_SUBSTATES_NUMBER;
  int ssNum_i = INT_SUBSTATES_NUMBER;
  size_t elNum = model->columns * model->rows * model->slices;

  long int outIndex = 0;
  int i;
  unsigned int j;

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    for (i = 0; i < ssNum_r; i++) {
      for (j = 0; j < elNum; j++)
        model->pQr_array[slot*ssNum_r+i]->current[j] = realSubstate_current_OUT[outIndex++];
    }

  outIndex = 0;

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    for (i = 0; i < ssNum_i; i++) {
      for (j = 0; j < elNum; j++)
        model->pQi_array[slot*ssNum_i+i]->current[j] = intSubstate_current_OUT[outIndex++];
    }
}

void realSubstatesMapper3D(struct CALModel3D *model, CALreal * current, CALreal * next) {
  //int ssNum = model->sizeof_pQr_array;
  int ssNum = REAL_SUBSTATES_NUMBER;
  size_t elNum = model->columns * model->rows * model->slices;
  long int outIndex = 0;
  long int outIndex1 = 0;
  int i;
  unsigned int j;

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    for (i = 0; i < ssNum; i++) {
        for (j = 0; j < elNum; j++)
          current[outIndex++] = model->pQr_array[slot*ssNum+i]->current[j];
        for (j = 0; j < elNum; j++)
          next[outIndex1++] = model->pQr_array[slot*ssNum+i]->next[j];
      }
}

void intSubstatesMapper3D(struct CALModel3D *model, CALint * current, CALint * next) {
  //int ssNum = model->sizeof_pQi_array;
  int ssNum = INT_SUBSTATES_NUMBER;
  size_t elNum = model->columns * model->rows * model->slices;
  long int outIndex = 0;\
  long int outIndex1 = 0;
  int i;
  unsigned int j;

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    for (i = 0; i < ssNum; i++) {
        for (j = 0; j < elNum; j++)
          current[outIndex++] = model->pQi_array[slot*ssNum+i]->current[j];
        for (j = 0; j < elNum; j++)
          next[outIndex1++] = model->pQi_array[slot*ssNum+i]->next[j];
      }
}


//------------------------------------------------------------------

void updateF(CALreal* Q_current,CALreal* Q_next)
{

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    {
      //calCopyBuffer3Dr(Q_Fx_next[slot],Q_Fx_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
      //calCopyBuffer3Dr(Q_Fy_next[slot],Q_Fy_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
      //calCopyBuffer3Dr(Q_Fz_next[slot],Q_Fz_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
//      #pragma omp target
//      calCopyBuffer3Dr(REAL_SUBSTATE(Q_next,slot,FX),REAL_SUBSTATE(Q_current,slot,FX), X_CELLS, Y_CELLS, Z_CELLS);
//#pragma omp target
//calCopyBuffer3Dr(REAL_SUBSTATE(Q_next,slot,FY),REAL_SUBSTATE(Q_current,slot,FY), X_CELLS, Y_CELLS, Z_CELLS);
//#pragma omp target
//calCopyBuffer3Dr(REAL_SUBSTATE(Q_next,slot,FZ),REAL_SUBSTATE(Q_current,slot,FZ), X_CELLS, Y_CELLS, Z_CELLS);
      CALint size=X_CELLS*Y_CELLS*Z_CELLS;
#ifdef GPU
#pragma omp target teams distribute parallel for num_teams(num_Teams) num_threads(num_Threads) dist_schedule(static,64) schedule(static,1)
#endif
   for(int i=0;i<size;i++)
    {
       REAL_SUBSTATE(Q_current,slot,FX)[i]=REAL_SUBSTATE(Q_next,slot,FX)[i];
       REAL_SUBSTATE(Q_current,slot,FY)[i]=REAL_SUBSTATE(Q_next,slot,FY)[i];
       REAL_SUBSTATE(Q_current,slot,FZ)[i]=REAL_SUBSTATE(Q_next,slot,FZ)[i];
    }
    }
}

void updateP(CALreal* Q_current,CALreal* Q_next)
{

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    {
      //calCopyBuffer3Dr(Q_px_next[slot],Q_px_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
      //calCopyBuffer3Dr(Q_py_next[slot],Q_py_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
      //calCopyBuffer3Dr(Q_pz_next[slot],Q_pz_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
//      #pragma omp target
//      calCopyBuffer3Dr(REAL_SUBSTATE(Q_next,slot,PX),REAL_SUBSTATE(Q_current,slot,PX), X_CELLS, Y_CELLS, Z_CELLS);
//      #pragma omp target
//      calCopyBuffer3Dr(REAL_SUBSTATE(Q_next,slot,PY),REAL_SUBSTATE(Q_current,slot,PY), X_CELLS, Y_CELLS, Z_CELLS);
//      #pragma omp target
//      calCopyBuffer3Dr(REAL_SUBSTATE(Q_next,slot,PZ),REAL_SUBSTATE(Q_current,slot,PZ), X_CELLS, Y_CELLS, Z_CELLS);
      CALint size=X_CELLS*Y_CELLS*Z_CELLS;
#ifdef GPU
#pragma omp target teams distribute parallel for num_teams(num_Teams) num_threads(num_Threads) dist_schedule(static,64) schedule(static,1)
#endif
   for(int i=0;i<size;i++)
    {
       REAL_SUBSTATE(Q_current,slot,PX)[i]=REAL_SUBSTATE(Q_next,slot,PX)[i];
       REAL_SUBSTATE(Q_current,slot,PY)[i]=REAL_SUBSTATE(Q_next,slot,PY)[i];
       REAL_SUBSTATE(Q_current,slot,PZ)[i]=REAL_SUBSTATE(Q_next,slot,PZ)[i];
    }


  }
}

void updateV(CALreal* Q_current,CALreal* Q_next)
{

    for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    {
      //calCopyBuffer3Dr(Q_vx_next[slot],Q_vx_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
      //calCopyBuffer3Dr(Q_vy_next[slot],Q_vy_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
      //calCopyBuffer3Dr(Q_vz_next[slot],Q_vz_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
//      #pragma omp target
//      calCopyBuffer3Dr(REAL_SUBSTATE(Q_next,slot,VX),REAL_SUBSTATE(Q_current,slot,VX), X_CELLS, Y_CELLS, Z_CELLS);
//      #pragma omp target
//      calCopyBuffer3Dr(REAL_SUBSTATE(Q_next,slot,VY),REAL_SUBSTATE(Q_current,slot,VY), X_CELLS, Y_CELLS, Z_CELLS);
//      #pragma omp target
//      calCopyBuffer3Dr(REAL_SUBSTATE(Q_next,slot,VZ),REAL_SUBSTATE(Q_current,slot,VZ), X_CELLS, Y_CELLS, Z_CELLS);
        CALint size=X_CELLS*Y_CELLS*Z_CELLS;
  #ifdef GPU
  #pragma omp target teams distribute parallel for num_teams(num_Teams) num_threads(num_Threads) dist_schedule(static,64) schedule(static,1)
  #endif
     for(int i=0;i<size;i++)
      {
         REAL_SUBSTATE(Q_current,slot,VX)[i]=REAL_SUBSTATE(Q_next,slot,VX)[i];
         REAL_SUBSTATE(Q_current,slot,VY)[i]=REAL_SUBSTATE(Q_next,slot,VY)[i];
         REAL_SUBSTATE(Q_current,slot,VZ)[i]=REAL_SUBSTATE(Q_next,slot,VZ)[i];
      }


    }
}

void updateID(CALint* ID_current,CALint* ID_next)
{

  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
  {
      //calCopyBuffer3Di(Q_ID_next[slot],Q_ID_current[slot], X_CELLS, Y_CELLS, Z_CELLS);
//      #pragma omp target
//    calCopyBuffer3Di(INT_SUBSTATE(ID_next,slot,PID),INT_SUBSTATE(ID_current,slot,PID), X_CELLS, Y_CELLS, Z_CELLS);

      CALint size=X_CELLS*Y_CELLS*Z_CELLS;
#ifdef GPU
#pragma omp target teams distribute parallel for num_teams(num_Teams) num_threads(num_Threads) dist_schedule(static,64) schedule(static,1)
#endif
   for(int i=0;i<size;i++)
      INT_SUBSTATE(ID_current,slot,PID)[i]=INT_SUBSTATE(ID_next,slot,PID)[i];

  }
}

void transferDataToDevice()
{
#pragma omp target enter data map(to:ID_current[0:ID_SIZE])
#pragma omp target enter data map(to:ID_next[0:ID_SIZE])
#pragma omp target enter data map(to:Q_current[0:Q_SIZE])
#pragma omp target enter data map(to:Q_next[0:Q_SIZE])
#pragma omp target enter data map(to:Xi[0:SIZE_OF_X])
#pragma omp target enter data map(to:Xj[0:SIZE_OF_X])
#pragma omp target enter data map(to:Xk[0:SIZE_OF_X])

}

void requestUpdateFromDevice()
{
#pragma omp target update from(ID_current[0:ID_SIZE])
#pragma omp target update from(ID_next[0:ID_SIZE])
#pragma omp target update from(Q_current[0:Q_SIZE])
#pragma omp target update from(Q_next[0:Q_SIZE])
//#pragma omp target exit data map(from:ID_current[0:ID_SIZE])
//#pragma omp target exit data map(from:ID_next[0:ID_SIZE])
//#pragma omp target exit data map(from:Q_current[0:Q_SIZE])
//#pragma omp target exit data map(from:Q_next[0:Q_SIZE])
}

void resetForce(CALint * ID_current,CALreal* Q_next)
{
#ifdef GPU
#pragma omp target teams distribute parallel for num_teams(num_Teams) num_threads(num_Threads) dist_schedule(static,64) schedule(static,1) collapse(3)
#endif
  for (int i = 0; i < X_CELLS; i++)
    for (int j = 0; j < Y_CELLS; j++)
      for (int k = 0; k < Z_CELLS; k++)
        resetF(ID_current,Q_next,i, j, k);

}

void innerCollision(CALint * ID_current,CALreal*Q_current,CALreal* Q_next)
{

#ifdef GPU
#pragma omp target teams distribute parallel for num_teams(num_Teams) num_threads(num_Threads) dist_schedule(static,64) schedule(static,1) collapse(3)
#endif
  for (int i = 0; i < X_CELLS; i++)
    for (int j = 0; j < Y_CELLS; j++)
      for (int k = 0; k < Z_CELLS; k++)
        inner_collision(ID_current,Q_current,Q_next, i, j, k);
}

void outerCollision(CALint * ID_current,CALreal*Q_current,CALreal* Q_next, CALint* Xi, CALint* Xj, CALint* Xk)
{

#ifdef GPU
#pragma omp target teams distribute parallel for num_teams(num_Teams) num_threads(num_Threads) dist_schedule(static,64) schedule(static,1) collapse(3)
#endif
  for (int i = 0; i < X_CELLS; i++)
    for (int j = 0; j < Y_CELLS; j++)
      for (int k = 0; k < Z_CELLS; k++)
        outer_collision(ID_current,Q_current,Q_next,Xi,Xj,Xk, i, j, k);
}

void muovili(CALint * ID_current,CALreal*Q_current,CALreal* Q_next)
{
#ifdef GPU
#pragma omp target teams distribute parallel for num_teams(num_Teams) num_threads(num_Threads) dist_schedule(static,64) schedule(static,1) collapse(3)
#endif
   for (int i = 0; i < X_CELLS; i++)
    for (int j = 0; j < Y_CELLS; j++)
      for (int k = 0; k < Z_CELLS; k++)
        movili(ID_current,Q_current,Q_next, i, j, k);
}

void muoviliCazzo(CALint * ID_current,CALint *ID_next,CALreal* Q_current,CALreal* Q_next, CALint* Xi, CALint* Xj, CALint* Xk)
{
#ifdef GPU
#pragma omp target teams distribute parallel for num_teams(num_Teams) num_threads(num_Threads) dist_schedule(static,64) schedule(static,1) collapse(3)
#endif
    for (int i = 0; i < X_CELLS; i++)
    for (int j = 0; j < Y_CELLS; j++)
      for (int k = 0; k < Z_CELLS; k++)
        moviliCazzu(ID_current,ID_next,Q_current,Q_next, Xi,Xj,Xk, i, j, k);

}

void transizioniGlobali(struct CALModel3D* modello)
{


resetForce(ID_current,Q_next);

//update
//#pragma omp target
updateF(Q_current,Q_next);


innerCollision(ID_current,Q_current,Q_next);
outerCollision(ID_current,Q_current,Q_next,Xi,Xj,Xk);


//update
//#pragma omp target
updateF(Q_current,Q_next);


muovili(ID_current,Q_current,Q_next);


//update
//#pragma omp target
updateP(Q_current,Q_next);
//#pragma omp target
updateV(Q_current,Q_next);


muoviliCazzo(ID_current,ID_next,Q_current,Q_next,Xi,Xj,Xk);

//update
//#pragma omp target
updateF(Q_current,Q_next);
//#pragma omp target
updateP(Q_current,Q_next);
//#pragma omp target
updateV(Q_current,Q_next);
//#pragma omp target
updateID(ID_current,ID_next);

  elapsed_time += DELTA_T;
#ifdef GPU
  if(step==STEPS)
  requestUpdateFromDevice();
#endif

#ifdef VERBOSE
  printSummary(modello);
#endif

  CALint S = INTEGRITY_CHECK_STEPS;
  if (step % S == 0)
    {
      CALint missing_particle = findMissingParticle(modello);
      if (missing_particle)
        {
          printf("ERROR: missing particle with ID %d\n", missing_particle);
          exit(EXIT_FAILURE);
        }
    }

}

CALbyte runCAStep3D(struct CALModel3D* modello)
{
  transizioniGlobali(modello);

  if (caminalu(modello))
    return CAL_FALSE;
  return CAL_TRUE;
}

void partilu()
{
  u_modellu = calCADef3D(X_CELLS,Y_CELLS,Z_CELLS,CAL_MOORE_NEIGHBORHOOD_3D,CAL_SPACE_TOROIDAL,CAL_NO_OPT);

  Q.Fx = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.Fy = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.Fz = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.px = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.py = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.pz = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.vx = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.vy = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.vz = (struct CALSubstate3Dr**)malloc(sizeof(struct CALSubstate3Dr*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q.ID = (struct CALSubstate3Di**)malloc(sizeof(struct CALSubstate3Di*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);

  for(int slot=0;slot<MAX_NUMBER_OF_PARTICLES_PER_CELL;slot++)
    {
      Q.Fx[slot] = calAddSubstate3Dr(u_modellu);
      Q.Fy[slot] = calAddSubstate3Dr(u_modellu);
      Q.Fz[slot] = calAddSubstate3Dr(u_modellu);
      Q.px[slot] = calAddSubstate3Dr(u_modellu);
      Q.py[slot] = calAddSubstate3Dr(u_modellu);
      Q.pz[slot] = calAddSubstate3Dr(u_modellu);
      Q.vx[slot] = calAddSubstate3Dr(u_modellu);
      Q.vy[slot] = calAddSubstate3Dr(u_modellu);
      Q.vz[slot] = calAddSubstate3Dr(u_modellu);
      Q.ID[slot] = calAddSubstate3Di(u_modellu);

      calInitSubstate3Dr(u_modellu,Q.Fx[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.Fy[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.Fz[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.px[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.py[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.pz[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.vx[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.vy[slot],0.0);
      calInitSubstate3Dr(u_modellu,Q.vz[slot],0.0);
      calInitSubstate3Di(u_modellu,Q.ID[slot],NULL_ID);
    }

/*
  Q_Fx_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_Fx_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_Fy_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_Fy_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_Fz_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_Fz_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_px_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_px_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_py_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_py_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_pz_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_pz_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_vx_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_vx_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_vy_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_vy_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_vz_current = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_vz_next    = (CALreal**)malloc(sizeof(CALreal*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_ID_current = (CALint**)malloc(sizeof(CALint*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);
  Q_ID_next    = (CALint**)malloc(sizeof(CALint*)*MAX_NUMBER_OF_PARTICLES_PER_CELL);

  for(int slot=0;slot<MAX_NUMBER_OF_PARTICLES_PER_CELL;slot++)
    {
      Q_Fx_current[slot] = Q.Fx[slot]->current;
      Q_Fx_next[slot]    = Q.Fx[slot]->next;
      Q_Fy_current[slot] = Q.Fy[slot]->current;
      Q_Fy_next[slot]    = Q.Fy[slot]->next;
      Q_Fz_current[slot] = Q.Fz[slot]->current;
      Q_Fz_next[slot]    = Q.Fz[slot]->next;
      Q_px_current[slot] = Q.px[slot]->current;
      Q_px_next[slot]    = Q.px[slot]->next;
      Q_py_current[slot] = Q.py[slot]->current;
      Q_py_next[slot]    = Q.py[slot]->next;
      Q_pz_current[slot] = Q.pz[slot]->current;
      Q_pz_next[slot]    = Q.pz[slot]->next;
      Q_vx_current[slot] = Q.vx[slot]->current;
      Q_vx_next[slot]    = Q.vx[slot]->next;
      Q_vy_current[slot] = Q.vy[slot]->current;
      Q_vy_next[slot]    = Q.vy[slot]->next;
      Q_vz_current[slot] = Q.vz[slot]->current;
      Q_vz_next[slot]    = Q.vz[slot]->next;
      Q_ID_current[slot] = Q.ID[slot]->current;
      Q_ID_next[slot]    = Q.ID[slot]->next;
    }
*/

  // Boundary
  boundaryCellsSerial(u_modellu);

  // Initial conditions
  initial_nummber_of_particles = 0;
  elapsed_time = 0.0;
  mmiscali_nta_cella_seriale(u_modellu);
  cancella_particelle_in_urto(u_modellu);

  // Simulation step
  step = 1;



  // Linearized substates
  Q_current  = (CALreal*)malloc(sizeof(CALreal)*Q_SIZE);
  Q_next     = (CALreal*)malloc(sizeof(CALreal)*Q_SIZE);
  ID_current = (CALint*)malloc(sizeof(CALint)  *ID_SIZE);
  ID_next    = (CALint*)malloc(sizeof(CALint)  *ID_SIZE);

  realSubstatesMapper3D(u_modellu, Q_current, Q_next);
  intSubstatesMapper3D(u_modellu, ID_current, ID_next);

  // allocate and build neighborhood arrays
  Xi = (CALint*)malloc(sizeof(CALint)*u_modellu->sizeof_X);
  Xj = (CALint*)malloc(sizeof(CALint)*u_modellu->sizeof_X);
  Xk = (CALint*)malloc(sizeof(CALint)*u_modellu->sizeof_X);
  for (int n=0; n<u_modellu->sizeof_X; n++)
    {
      Xi[n] = u_modellu->X[n].i;
      Xj[n] = u_modellu->X[n].j;
      Xk[n] = u_modellu->X[n].k;
    }
#ifdef GPU
    transferDataToDevice();
#endif


#ifdef VERBOSE
  printf("The 3D particles computational model\n");
#ifdef OMP
  printf("OpenMP parallel execution enabled!\n");
#endif
  printf("A system of %d particles will be simulated for %f s, subdivided in %d steps, each one corresponding to %f s\n", initial_nummber_of_particles, TOTAL_SIMULATION_TIME, STEPS, DELTA_T);
#endif
}
