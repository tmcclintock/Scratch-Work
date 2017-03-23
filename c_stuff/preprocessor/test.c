#include <stdio.h>
#include <time.h>

int main(){
  int x = 0;
  time_t rawtime = time(&rawtime);

  printf("real print statement; x=%d\n",x);
  //printf("time = %s\n",ctime(&rawtime));

  #ifdef TEST
  x+=1;
  printf("real print statement2; x=%d\n",x);
  #endif

  #ifdef LOL
  x+=1;
  printf("real print statement3; x=%d\n",x);
  #endif

}
