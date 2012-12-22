#include "mex.h" /* Always include this */

void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
int nrhs, const mxArray *prhs[]) /* Input variables */
{
    int i;
mexPrintf("Salam zendegi\n"); /* Do something interesting */
mexPrintf("%d\n",i);
}