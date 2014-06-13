#include <stdlib.h>

int SYSTEM(const char *command )
{
  return system (command ); 
}

void ABORT(void)
{
  abort();
}
