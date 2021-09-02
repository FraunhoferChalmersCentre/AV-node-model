#include "mex.h"
#include <math.h>
#include <vector>
#include <queue>

using namespace std;

inline double res(double DI,double *p)  { return p[0]+p[1]*(1.0-exp(-DI/p[2])); }
inline double res2(double DI,double *p) { return p[0]+p[1]*exp(-DI/p[2]); }

// [out, pos_AT, AT_out, everything] = run_new_model2(impulse,apd_f,apd_s,apd_n,delay_f,delay_s,delay_n) 

struct imp
{
    double time;
    int ind;
    int impind;
    int pathOg;
    int pathEnd;
};

class cmpimp {
public:
    bool operator()(imp& i1, imp& i2)
    {
        if (i1.time > i2.time) return true;
        return false;
    }
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *impulse =   (double *) mxGetPr(prhs[0]);
    double *apd_f =     (double *) mxGetPr(prhs[1]);
    double *apd_s =     (double *) mxGetPr(prhs[2]);
    double  apd_n_s =              mxGetScalar(prhs[3]);
    double *delay_f =   (double *) mxGetPr(prhs[4]);
    double *delay_s =   (double *) mxGetPr(prhs[5]);
    double *RT_in =     (double *) mxGetPr(prhs[6]); 
    int     MaxLoopsInt =          mxGetScalar(prhs[7]);
    
    int nRT_in = mxGetNumberOfElements(prhs[6]);
    int nImpulses = mxGetNumberOfElements(prhs[0]);
    
    
    double* apd_n = (double *) mxCalloc(3, sizeof(double));
    apd_n[0] = apd_n_s;
    apd_n[1] = 0;
    apd_n[2] = 0;
    
    plhs[0] = mxCreateDoubleMatrix(nImpulses*2,1,mxREAL);
    plhs[12] = mxCreateDoubleMatrix(21,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nImpulses*2,1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nImpulses*2,1,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(21,1,mxREAL);
    plhs[5] = mxCreateDoubleMatrix(nImpulses*2*2,1,mxREAL);
    plhs[6] = mxCreateDoubleMatrix(nImpulses*2*2,1,mxREAL);
    plhs[7] = mxCreateDoubleMatrix(nImpulses*2*2,1,mxREAL);

    plhs[8] = mxCreateDoubleMatrix(nImpulses*2*21,1,mxREAL);
    plhs[9] = mxCreateDoubleMatrix(nImpulses*2*21,1,mxREAL);
    plhs[10] = mxCreateDoubleMatrix(nImpulses*2*21,1,mxREAL);
    plhs[11] = mxCreateDoubleMatrix(nImpulses*2*21,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nImpulses*2*21,1,mxREAL);
    
    double* out = mxGetPr(plhs[0]);
    
    double* pathwayEnd = mxGetPr(plhs[2]);
    double* pathway = mxGetPr(plhs[3]);    
    double* ATnode = mxGetPr(plhs[5]);
    double* ATtime = mxGetPr(plhs[6]);
    double* ATpath = mxGetPr(plhs[7]);
    
    double* ATnode2 = mxGetPr(plhs[8]);
    double* ATtime2 = mxGetPr(plhs[9]);
    double* ATpath2 = mxGetPr(plhs[10]);
    double* ATimp2 = mxGetPr(plhs[11]);
    
    double* ListOfTime = mxGetPr(plhs[1]);
    
    double* DL = mxGetPr(plhs[4]);
    double* AT = (double *) mxCalloc(21, sizeof(double));
    
    double* RT = mxGetPr(plhs[12]);
    
    for (int i=0; i<nRT_in; i++) RT[i] = RT_in[i];
    
    int current_index(0);
    
    int i(0), c(0), counter(0), d(0), e(0), numberOfActivations(0);
    int next_pos, origin, pathOg, pathEnd;
    long int g(0);
    double next_time;
    
    priority_queue<imp, vector<imp>, cmpimp> q;
    
    imp s;
    for (int i = 0; i < nImpulses; ++i)
    {
        s.time = impulse[i];
        s.impind = i+1;
        s.ind = 0;
        s.pathOg = 1;
        s.pathEnd = 1;
        q.push(s);
        s.ind = 10;
        s.pathOg = 2;
        s.pathEnd = 2;
        q.push(s);
    }
    
    while (!q.empty())
    {
        if ( c > (MaxLoopsInt) )
        {
            break;
        }
        
        s = q.top();
        q.pop();
                
        next_time = s.time;
        next_pos = s.ind;
        origin = s.impind;   
        pathOg = s.pathOg;
        pathEnd = s.pathEnd;

        ListOfTime[g] = next_time-RT[next_pos];
        g = g + 1;
         
        // check delay thresholds
        double tmp; 
        bool flag{true};
        
        if (res(next_time-RT[next_pos],apd_s)<0) flag = false;
        
        if ( (next_time >= RT[next_pos]) & flag)
        {
            // If it gets stuck
            numberOfActivations = numberOfActivations + 1;
            if (numberOfActivations > MaxLoopsInt*21*20)
            {
                break;
            }
            //
            if (next_pos<10) // slow pathway
            {
      
                DL[next_pos] = res2(next_time-RT[next_pos],delay_s);
                RT[next_pos] = next_time + res(next_time-RT[next_pos],apd_s);
                AT[next_pos] = next_time;
                
                
                
                if (next_pos==9)
                {
                    s.time = AT[next_pos]+DL[next_pos];
                    s.impind = origin;
                    s.ind = 8;
                    q.push(s);
                    s.ind = 19;
                    q.push(s); 
                    s.ind = 20;
                    s.pathEnd = 1;
                    q.push(s);
                }
                else
                {
                    s.time = AT[next_pos]+DL[next_pos];
                    s.impind = origin;
                    s.ind = next_pos+1;
                    q.push(s);
                    if (next_pos>0)
                    {
                        s.ind = next_pos-1;
                        q.push(s);
                    }
                }
            }
            else if (next_pos<20) // fast pathway
            {
                DL[next_pos] = res2(next_time-RT[next_pos],delay_f);
                RT[next_pos] = next_time + res(next_time-RT[next_pos],apd_f);
                AT[next_pos] = next_time;

                if (next_pos==19)
                {
                    s.time = AT[next_pos]+DL[next_pos];
                    s.impind = origin;
                    s.ind = 18;
                    q.push(s);
                    s.ind = 9;
                    q.push(s);
                    s.ind = 20;
                    s.pathEnd = 2;
                    q.push(s);
                }
                else
                {
                    s.time = AT[next_pos]+DL[next_pos];
                    s.impind = origin;
                    s.ind = next_pos+1;
                    q.push(s);
                    if (next_pos>10)
                    {
                        s.ind = next_pos-1;
                        q.push(s);
                    }
                }
            }
            else if (next_pos==20) // nodal-his 
            {
                RT[next_pos] = next_time + res(next_time-RT[next_pos],apd_n);
                AT[next_pos] = next_time;   
                if (c==2*nImpulses) break;
                out[c] = AT[next_pos]+60;
                pathway[c] = pathOg;
                pathwayEnd[c] = pathEnd;
                c = c + 1;
            }
            
            ATnode2[e] = next_pos;
            ATtime2[e] = next_time;
            ATpath2[e] = pathOg;
            ATimp2[e] = origin;
            e++;
                   
        } // end if
        else if (d<nImpulses*2*2) 
        {
            ATnode[d] = next_pos;
            ATtime[d] = next_time;
            ATpath[d] = pathOg;
            d++;
        }
        
        
    } // end while
    mxFree(AT);
    mxFree(apd_n);
}


