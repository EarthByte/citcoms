#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

int Emergency_stop;

void parallel_process_termination();



void interuption()
{
  if (Emergency_stop)
    parallel_process_termination();
  else
    Emergency_stop++;
  fprintf(stderr,"Cleaning up before exit\n");
  return;
}


void set_signal()
{
  Emergency_stop = 0;

  signal(SIGINT,interuption);
  signal(SIGTERM,interuption);
  return;
}



