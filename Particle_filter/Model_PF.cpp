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
    double *apd_n =     (double *) mxGetPr(prhs[3]);
    double *delay_f =   (double *) mxGetPr(prhs[4]);
    double *delay_s =   (double *) mxGetPr(prhs[5]);
    double *delay_n =   (double *) mxGetPr(prhs[6]);
    double *RT_in =     (double *) mxGetPr(prhs[7]);
    double *q_in =     (double *) mxGetPr(prhs[8]);
    
    int nRT_in = mxGetNumberOfElements(prhs[7]);
    int nImpulses = mxGetNumberOfElements(prhs[0]);
    int nSave = (nImpulses + 100)*20+5000;
    int n_q_in = mxGetNumberOfElements(prhs[8]);
    
    plhs[0] = mxCreateDoubleMatrix(nSave,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nSave,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nSave,1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nSave,1,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(21,1,mxREAL);
    plhs[5] = mxCreateDoubleMatrix(nSave,1,mxREAL);
    plhs[6] = mxCreateDoubleMatrix(nSave,1,mxREAL);
    plhs[7] = mxCreateDoubleMatrix(nSave,1,mxREAL);

    plhs[8] = mxCreateDoubleMatrix(nSave,1,mxREAL);
    plhs[9] = mxCreateDoubleMatrix(nSave,1,mxREAL);
    plhs[10] = mxCreateDoubleMatrix(nSave,1,mxREAL);
    plhs[11] = mxCreateDoubleMatrix(nSave,1,mxREAL);
    plhs[12] = mxCreateDoubleMatrix(21,1,mxREAL);
    plhs[13] = mxCreateDoubleMatrix(nSave,1,mxREAL);
    plhs[14] = mxCreateDoubleMatrix(nSave,1,mxREAL);
    plhs[15] = mxCreateDoubleMatrix(nSave,1,mxREAL);
    plhs[16] = mxCreateDoubleMatrix(5, 10000, mxREAL);
    plhs[17] = mxCreateDoubleMatrix(2, nSave, mxREAL);
    plhs[18] = mxCreateDoubleMatrix(2, nSave, mxREAL);
    
    
    //double *DI = mxGetPr(plhs[4]);
    double* out = mxGetPr(plhs[0]);
    double* ListOfTime = mxGetPr(plhs[1]);
    double* pathwayEnd = mxGetPr(plhs[2]);
    
    double* pathway = mxGetPr(plhs[3]);
        
    double* ATnode = mxGetPr(plhs[5]);
    double* ATtime = mxGetPr(plhs[6]);
    double* ATpath = mxGetPr(plhs[7]);
    
    double* ATnode2 = mxGetPr(plhs[8]);
    double* ATtime2 = mxGetPr(plhs[9]);
    double* ATpath2 = mxGetPr(plhs[10]);
    double* ATimp2 = mxGetPr(plhs[11]);
    
    double* ListOfPos = mxGetPr(plhs[13]);
    double* ListOfD = mxGetPr(plhs[14]);
    double* ListOfR = mxGetPr(plhs[15]);

    double* q_out = mxGetPr(plhs[16]);

    double* FP_out = mxGetPr(plhs[17]);
    double* SP_out = mxGetPr(plhs[18]);
    
    //double* DL = (double *) mxCalloc(21, sizeof(double));
    double* DL = mxGetPr(plhs[4]);
    double* AT = (double *) mxCalloc(21, sizeof(double));
    
    double* RT = mxGetPr(plhs[12]);
    
    for (int i=0; i<nRT_in; i++) RT[i] = RT_in[i];
    
    //int* index = (int *) mxCalloc(21, sizeof(int));
    int current_index(0);
    
    int i(0), c(0), counter(0), d(0), e(0), numberOfActivations(0), c_queue(0), c_SP(0), c_FP(0);
    int next_pos, origin, pathOg, pathEnd;
    long int g(0);
    double next_time;
    
    priority_queue<imp, vector<imp>, cmpimp> q;

    imp s;
    for (int i = 0; i < nImpulses; ++i)
    {
        s.time = impulse[i];
        s.impind = 0;
        s.ind = 0;
        s.pathOg = 1;
        s.pathEnd = 1;
        q.push(s);
        //s.time = impulse[i]+20;
        s.ind = 10;
        s.pathOg = 2;
        s.pathEnd = 2;
        q.push(s);
    }
    
    for (int i_q_in = 0; i_q_in < n_q_in/5; ++i_q_in)
    {
        s.time      = q_in[i_q_in*5];
        s.ind       = q_in[i_q_in*5 + 1];
        s.impind    = q_in[i_q_in*5 + 2];
        s.pathOg    = q_in[i_q_in*5 + 3];
        s.pathEnd   = q_in[i_q_in*5 + 4];
        //mexPrintf("%f ", q_in[i_q_in]);
        q.push(s);
    }


    
    while (!q.empty())
    {        
        s = q.top();
        q.pop();
                
        next_time = s.time;
        next_pos = s.ind;
        //origin = s.impind;   
        pathOg = s.pathOg;
        pathEnd = s.pathEnd;

        ListOfTime[g] = next_time;
        ListOfPos[g] = next_pos + 1;
        if (next_pos<10)
        {
            ListOfR[g] = res(next_time-RT[next_pos],apd_s);
            ListOfD[g] = res2(next_time-RT[next_pos],delay_s);
        }
        else if (next_pos<20)
        {
            ListOfR[g] = res(next_time-RT[next_pos],apd_f);
            ListOfD[g] = res2(next_time-RT[next_pos],delay_f);
        }
        else if (next_pos == 20)
        {
            ListOfR[g] = res(next_time-RT[next_pos],apd_n);
            ListOfD[g] = res2(next_time-RT[next_pos],delay_n);
        }
        g = g + 1;
         
        double tmp; 
        
        if (next_time >= RT[next_pos])
        {
            // If it get stuck
            numberOfActivations = numberOfActivations + 1;
            if (numberOfActivations > 1000)
            {
                break;
            }
            //
            if (next_pos<10) // slow pathway
            {
      
                DL[next_pos] = res2(next_time-RT[next_pos],delay_s);
                RT[next_pos] = next_time + res(next_time-RT[next_pos],apd_s);
                AT[next_pos] = next_time;

                SP_out[2*c_SP] =  DL[next_pos];
                SP_out[2*c_SP + 1] =  RT[next_pos] - next_time;
                c_SP = c_SP + 1;

                if (next_pos==9)
                {
                    s.time = AT[next_pos]+DL[next_pos];
                    s.impind = 1;
                    s.ind = 8;
                    q.push(s);
                    s.ind = 19;
                    q.push(s); 
                    s.ind = 20;
                    s.pathEnd = 1;
                    q.push(s);
                    //mexPrintf("push from slow\n");
                    //out[c] = AT[next_pos];
                    //c = c + 1;
                }
                else
                {
                    s.time = AT[next_pos]+DL[next_pos];
                    s.impind = 1;
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

                FP_out[2*c_FP] =  DL[next_pos];
                FP_out[2*c_FP + 1] =  RT[next_pos] - next_time;
                c_FP = c_FP + 1;

                //q = [q AT(next_pos)+DL(next_pos) AT(next_pos)+DL(next_pos)];
                if (next_pos==19)
                {
                    s.time = AT[next_pos]+DL[next_pos];
                    s.impind = 1;
                    s.ind = 18;
                    q.push(s);
                    s.ind = 9;
                    q.push(s);
                    s.ind = 20;
                    s.pathEnd = 2;
                    q.push(s);
                    //out[c] = AT[next_pos];
                    //c = c + 1;
                    //out(c,:) = [AT(next_pos) i];
                    //c = c + 1;
                    //mexPrintf("push from fast\n");
                }
                else
                {
                    s.time = AT[next_pos]+DL[next_pos];
                    s.impind = 1;
                    s.ind = next_pos+1;
                    q.push(s);
                    if (next_pos>10)
                    {
                        s.ind = next_pos-1;
                        q.push(s);
                    }
                }
            }
            else if (next_pos==20) // central node
            {
                DL[next_pos] = res2(next_time-RT[next_pos],delay_n);
                RT[next_pos] = next_time + res(next_time-RT[next_pos],apd_n);
                AT[next_pos] = next_time;   

                out[c] = AT[next_pos] + DL[next_pos];
                //out_ind[c] = origin;
                pathway[c] = pathOg;
                pathwayEnd[c] = pathEnd;
                c = c + 1;

                // run the queue out to save it
                while (!q.empty())
                {
                    s = q.top();
                    q.pop();
                            
                    q_out[c_queue*5] = s.time;
                    q_out[c_queue*5 + 1] = s.ind;
                    q_out[c_queue*5 + 2] = s.impind;   
                    q_out[c_queue*5 + 3] = s.pathOg;
                    q_out[c_queue*5 + 4] = s.pathEnd;

                    //q_out[1001] = 1001;

                    c_queue = c_queue + 1;
            
                }
                break;
            }
            
            ATnode2[e] = next_pos;
            ATtime2[e] = next_time;
            ATpath2[e] = pathOg;
            ATimp2[e] = origin;
            e++;
            
            //mexPrintf("%g %u\n",next.time,next.ind);
        
        } // end if
        else if (d<1000)
        {
            ATnode[d] = next_pos;
            ATtime[d] = next_time;
            ATpath[d] = pathOg;
            d++;
        }
        
        
    } // end while
    //mxFree(DL);
    mxFree(AT);

    //mxFree(index);
}


