#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main () {
   int val;
//    char ps[200];
    // char *str;
   char *ps = (char *)malloc(200*sizeof(char));

   strcpy(ps, "V     -1 0 0 0 0 0 2   290 0");
   const char *str = strtok(ps, " ");
//    printf("%s \n", strchr(str+1, ' '));
//    printf("============== \n");
//    str = strchr(str+1, ' ');
//    printf("%s \n", str);
//    str = strchr(str+1, ' ');
//    printf("%s \n", str);
//    str = strchr(str+1, ' ');
//    printf("%s \n", str);
//    str = strchr(str+1, ' ');
//    printf("%s \n", str);
//    str = strchr(str+1, ' ');
//    printf("%s \n", str);
//    str = strchr(str+1, ' ');
//    return(0);

    while (str != NULL) {
        printf("%s\n", str);
        str = strtok(NULL, " ");
    }

    return(0);
}