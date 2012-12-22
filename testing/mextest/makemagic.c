/* Run as "makemagic N" to make an NxN magic matrix */
#include <stdlib.h>
#include <stdio.h>
int main(int argc, char *argv[])
{
float v;
int n, N = (argc > 1) ? atoi(argv[1]) : 6;
/* Write a MATLAB script */
FILE *fp = fopen("makemagic.m", "wt");
fprintf(fp, "x = magic(%d);\n" /* Make an NxN magic matrix */
"save magic.txt x -ascii\n" /* Save to magic.txt */
"exit;", N); /* Exit MATLAB */
fclose(fp);
/* Call MATLAB to run the script */
printf("Calling MATLAB...\n");
system("matlab -r makemagic");
/* Read from the output file */
fp = fopen("magic.txt", "rt");
while(!feof(fp))
{
for(n = 0; n < N && fscanf(fp, "%f", &v) == 1; n++)
printf("%4d ", (int)v);
printf("\n");
}
fclose(fp);
return 0;
}