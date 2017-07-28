#include <ep_movili_cazzu.h>
#include <ep_utils.h>
#include <stdlib.h>

void pezziala(CALint* ID_next,CALreal* Q_next,int slot,int cell_x, int cell_y, int cell_z)
{
  //calSet3Dr(ca,Q.Fx[slot],cell_x,cell_y,cell_z,0.0);
  //calSet3Dr(ca,Q.Fy[slot],cell_x,cell_y,cell_z,0.0);
  //calSet3Dr(ca,Q.Fz[slot],cell_x,cell_y,cell_z,0.0);
  //calSet3Dr(ca,Q.px[slot],cell_x,cell_y,cell_z,0.0);
  //calSet3Dr(ca,Q.py[slot],cell_x,cell_y,cell_z,0.0);
  //calSet3Dr(ca,Q.pz[slot],cell_x,cell_y,cell_z,0.0);
  //calSet3Dr(ca,Q.vx[slot],cell_x,cell_y,cell_z,0.0);
  //calSet3Dr(ca,Q.vy[slot],cell_x,cell_y,cell_z,0.0);
  //calSet3Dr(ca,Q.vz[slot],cell_x,cell_y,cell_z,0.0);
  //calSet3Di(ca,Q.ID[slot],cell_x,cell_y,cell_z,NULL_ID);
  //calSetBuffer3DElement(Q_Fx_next[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  //calSetBuffer3DElement(Q_Fy_next[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  //calSetBuffer3DElement(Q_Fz_next[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  //calSetBuffer3DElement(Q_px_next[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  //calSetBuffer3DElement(Q_py_next[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  //calSetBuffer3DElement(Q_pz_next[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  //calSetBuffer3DElement(Q_vx_next[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  //calSetBuffer3DElement(Q_vy_next[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  //calSetBuffer3DElement(Q_vz_next[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  //calSetBuffer3DElement(Q_ID_next[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,NULL_ID);
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,slot,FX),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,slot,FY),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,slot,FZ),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,slot,PX),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,slot,PY),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,slot,PZ),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,slot,VX),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,slot,VY),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,slot,VZ),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,0.0);
  calSetBuffer3DElement(INT_SUBSTATE(ID_next,slot,PID),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,NULL_ID);
}

void sucala(CALint* ID_current,CALint* ID_next,CALreal* Q_current,CALreal* Q_next,CALint* Xi,CALint* Xj,CALint* Xk,int destination_slot, int source_slot, int cell_x, int cell_y, int cell_z, int n)
{
  //calSet3Dr(ca,Q.Fx[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.Fx[source_slot],cell_x,cell_y,cell_z,n));
  //calSet3Dr(ca,Q.Fy[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.Fy[source_slot],cell_x,cell_y,cell_z,n));
  //calSet3Dr(ca,Q.Fz[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.Fz[source_slot],cell_x,cell_y,cell_z,n));
  //calSet3Dr(ca,Q.px[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.px[source_slot],cell_x,cell_y,cell_z,n));
  //calSet3Dr(ca,Q.py[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.py[source_slot],cell_x,cell_y,cell_z,n));
  //calSet3Dr(ca,Q.pz[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.pz[source_slot],cell_x,cell_y,cell_z,n));
  //calSet3Dr(ca,Q.vx[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.vx[source_slot],cell_x,cell_y,cell_z,n));
  //calSet3Dr(ca,Q.vy[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.vy[source_slot],cell_x,cell_y,cell_z,n));
  //calSet3Dr(ca,Q.vz[destination_slot],cell_x,cell_y,cell_z, calGetX3Dr(ca,Q.vz[source_slot],cell_x,cell_y,cell_z,n));
  //calSet3Di(ca,Q.ID[destination_slot],cell_x,cell_y,cell_z, calGetX3Di(ca,Q.ID[source_slot],cell_x,cell_y,cell_z,n));
  //calSetBuffer3DElement(Q_Fx_next[destination_slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
  //   calGetBuffer3DElement(Q_Fx_current[source_slot],X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  //calSetBuffer3DElement(Q_Fy_next[destination_slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
  //   calGetBuffer3DElement(Q_Fy_current[source_slot],X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  //calSetBuffer3DElement(Q_Fz_next[destination_slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
  //   calGetBuffer3DElement(Q_Fz_current[source_slot],X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  //calSetBuffer3DElement(Q_px_next[destination_slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
  //   calGetBuffer3DElement(Q_px_current[source_slot],X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  //calSetBuffer3DElement(Q_py_next[destination_slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
  //   calGetBuffer3DElement(Q_py_current[source_slot],X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  //calSetBuffer3DElement(Q_pz_next[destination_slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
  //   calGetBuffer3DElement(Q_pz_current[source_slot],X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  //calSetBuffer3DElement(Q_vx_next[destination_slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
  //   calGetBuffer3DElement(Q_vx_current[source_slot],X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  //calSetBuffer3DElement(Q_vy_next[destination_slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
  //   calGetBuffer3DElement(Q_vy_current[source_slot],X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  //calSetBuffer3DElement(Q_vz_next[destination_slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
  //   calGetBuffer3DElement(Q_vz_current[source_slot],X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  //calSetBuffer3DElement(Q_ID_next[destination_slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
  //   calGetBuffer3DElement(Q_ID_current[source_slot],X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,destination_slot,FX),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
     calGetBuffer3DElement(REAL_SUBSTATE(Q_current,source_slot,FX),X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,destination_slot,FY),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
     calGetBuffer3DElement(REAL_SUBSTATE(Q_current,source_slot,FY),X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,destination_slot,FZ),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
     calGetBuffer3DElement(REAL_SUBSTATE(Q_current,source_slot,FZ),X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,destination_slot,PX),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
     calGetBuffer3DElement(REAL_SUBSTATE(Q_current,source_slot,PX),X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,destination_slot,PY),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
     calGetBuffer3DElement(REAL_SUBSTATE(Q_current,source_slot,PY),X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,destination_slot,PZ),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
     calGetBuffer3DElement(REAL_SUBSTATE(Q_current,source_slot,PZ),X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,destination_slot,VX),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
     calGetBuffer3DElement(REAL_SUBSTATE(Q_current,source_slot,VX),X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,destination_slot,VY),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
     calGetBuffer3DElement(REAL_SUBSTATE(Q_current,source_slot,VY),X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  calSetBuffer3DElement(REAL_SUBSTATE(Q_next,destination_slot,VZ),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
     calGetBuffer3DElement(REAL_SUBSTATE(Q_current,source_slot,VZ),X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
  calSetBuffer3DElement(INT_SUBSTATE(ID_next,destination_slot,PID),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z,
     calGetBuffer3DElement(INT_SUBSTATE(ID_current,source_slot,PID),X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)));
}

void moviliCazzu(CALint* ID_current,CALint* ID_next,CALreal* Q_current,CALreal* Q_next,CALint* Xi,CALint* Xj,CALint* Xk, int cell_x, int cell_y, int cell_z)
{
  CALreal px, py, pz;
  CALint  new_cell_x, new_cell_y, new_cell_z;

  //pezziali
  for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
    //if (calGet3Di(ca, Q.ID[slot],cell_x,cell_y,cell_z) > NULL_ID)
    //if (calGetBuffer3DElement(Q_ID_current[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z) > NULL_ID)
    if (calGetBuffer3DElement(INT_SUBSTATE(ID_current, slot, PID),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z) > NULL_ID)
      {
        //px = calGet3Dr(ca,Q.px[slot],cell_x,cell_y,cell_z);
        //py = calGet3Dr(ca,Q.py[slot],cell_x,cell_y,cell_z);
        //pz = calGet3Dr(ca,Q.pz[slot],cell_x,cell_y,cell_z);
        //px = calGetBuffer3DElement(Q_px_current[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z);
        //py = calGetBuffer3DElement(Q_py_current[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z);
        //pz = calGetBuffer3DElement(Q_pz_current[slot],X_CELLS,Y_CELLS,cell_x,cell_y,cell_z);
        px = calGetBuffer3DElement(REAL_SUBSTATE(Q_current,slot,PX),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z);
        py = calGetBuffer3DElement(REAL_SUBSTATE(Q_current,slot,PY),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z);
        pz = calGetBuffer3DElement(REAL_SUBSTATE(Q_current,slot,PZ),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z);

        new_cell_x = px/CELL_SIDE;
        new_cell_y = py/CELL_SIDE;
        new_cell_z = pz/CELL_SIDE;

        if ((cell_x != new_cell_x) || (cell_y != new_cell_y) || (cell_z != new_cell_z))
          pezziala(ID_next,Q_next,slot,cell_x,cell_y,cell_z);
      }

  //sucali
  for (int n=1; n<SIZE_OF_X; n++)
    for (int source_slot = 0; source_slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; source_slot++)
      //if (calGetX3Di(ca,Q.ID[source_slot],cell_x,cell_y,cell_z,n) > NULL_ID)
      //if (calGetBuffer3DElement(Q_ID_current[source_slot],X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)) > NULL_ID)
      if (calGetBuffer3DElement(INT_SUBSTATE(ID_current, source_slot, PID),X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS)) > NULL_ID)
        {
          //px = calGetX3Dr(ca,Q.px[source_slot],cell_x,cell_y,cell_z,n);
          //py = calGetX3Dr(ca,Q.py[source_slot],cell_x,cell_y,cell_z,n);
          //pz = calGetX3Dr(ca,Q.pz[source_slot],cell_x,cell_y,cell_z,n);
          //px = calGetBuffer3DElement(Q_px_current[source_slot],X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS));
          //py = calGetBuffer3DElement(Q_py_current[source_slot],X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS));
          //pz = calGetBuffer3DElement(Q_pz_current[source_slot],X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS));
          px = calGetBuffer3DElement(REAL_SUBSTATE(Q_current,source_slot,PX),X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS));
          py = calGetBuffer3DElement(REAL_SUBSTATE(Q_current,source_slot,PY),X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS));
          pz = calGetBuffer3DElement(REAL_SUBSTATE(Q_current,source_slot,PZ),X_CELLS,Y_CELLS,calGetToroidalX(cell_x+Xi[n],X_CELLS),calGetToroidalX(cell_y+Xj[n],Y_CELLS),calGetToroidalX(cell_z+Xk[n],Z_CELLS));

          new_cell_x = px/CELL_SIDE;
          new_cell_y = py/CELL_SIDE;
          new_cell_z = pz/CELL_SIDE;

          if ((cell_x == new_cell_x) && (cell_y == new_cell_y) && (cell_z == new_cell_z))
            {
              CALbyte sucked = CAL_FALSE;
              int destination_slot;
              for (destination_slot = 0; destination_slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; destination_slot++)
                //if (calGetNext3Di(ca,Q.ID[destination_slot],cell_x,cell_y,cell_z) == NULL_ID)
                if (calGetBuffer3DElement(INT_SUBSTATE(ID_next, destination_slot, PID),X_CELLS,Y_CELLS,cell_x,cell_y,cell_z) == NULL_ID)
                  {
                    sucala(ID_current,ID_next,Q_current,Q_next,Xi,Xj,Xk,destination_slot,source_slot,cell_x,cell_y,cell_z,n);
                    sucked = CAL_TRUE;
                    break;
                  }
              if (sucked == CAL_FALSE)
                {
                  printf("ERROR: unable to suck a particle.\n");
#ifdef VERBOSE
                  printf("cell_capacity: %d\n", MAX_NUMBER_OF_PARTICLES_PER_CELL);
                  printf("current_step: %d\n", step);
                  printf("source_slot: %d\n", source_slot);
                  printf("destination_slot: %d\n", destination_slot);
#endif
                 // exit(EXIT_FAILURE); problema nel device
                }
            }
        }
}
