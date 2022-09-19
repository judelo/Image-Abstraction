/*--------------------------- Commande MegaWave -----------------------------*/
/* mwcommand
  name = {sconst};
  version = {"1.0"};
  author = {"Jacques Froment"};
  function = {"Create a constant fsignal"};
  usage = {
  's':[size=512]->size    "number of samples in the signal (default 512)",
  'a':[amplitude=1.0]->A  "amplitude (default 1.0)",
  sconst<-signal          "output fsignal"
  };
*/
/*--- MegaWave - Copyright (C) 1992 Jacques Froment. All Rights Reserved. ---*/

#include <stdio.h>
#include <math.h>

/* Include always the MegaWave2 include file */
#include "mw3.h"

void sconst(int *size,
            float *A,
            Fsignal signal)

{

    signal = mw_change_fsignal(signal, *size);
    if (signal == NULL) mwerror(FATAL,1,"Not enough memory.");
    strcpy(signal->cmt,"Constant");

    mw_clear_fsignal(signal,*A);
}







