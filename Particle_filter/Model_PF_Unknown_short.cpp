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
    int AA_or_not;
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
    double *x_new =     (double *) mxGetPr(prhs[1]);
    double *x_old =     (double *) mxGetPr(prhs[2]);
    double *apd_n =     (double *) mxGetPr(prhs[3]);
    double *delay_n =   (double *) mxGetPr(prhs[4]);
    double *RT_in =     (double *) mxGetPr(prhs[5]);
    double *q_in =      (double *) mxGetPr(prhs[6]);
    int N_RR = static_cast<int>(mxGetScalar(prhs[7]));

    double *apd_f = new double[3];
    double *apd_s = new double[3];
    double *delay_f = new double[3];
    double *delay_s = new double[3];

    for (int i_3 = 0; i_3 < 3; ++i_3) {
    apd_f[i_3] = x_old[i_3];
    apd_s[i_3] = x_old[i_3+3];
    delay_f[i_3] = x_old[i_3+6];
    delay_s[i_3] = x_old[i_3+9];
    }

    int nImpulses = mxGetNumberOfElements(prhs[0]);
    int nSave = 5000;
    int n_q_in = mxGetNumberOfElements(prhs[7]);
    
    plhs[0] = mxCreateDoubleMatrix(nSave,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(2, nSave, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(2, nSave, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1, 2, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(21,nSave,mxREAL);



    double* out = mxGetPr(plhs[0]);
    double* FP_out = mxGetPr(plhs[1]);
    double* SP_out = mxGetPr(plhs[2]);
    double* AA_ind_out = mxGetPr(plhs[4]);
    double* RT_big = mxGetPr(plhs[5]);

    double DL[21];
    double* AT = (double *) mxCalloc(21, sizeof(double));
    double RT[21];
   
    
    for (int i=0; i<21; i++) RT[i] = RT_in[i];
       
    int i(0), c(0), numberOfActivations(0), c_SP(0), c_FP(0), i_numb(0);
    int next_pos, pathOg, pathEnd;
    double next_time;

    // Create a vector to store multiple versions of queues
    vector<priority_queue<imp, vector<imp>, cmpimp>> q_AA_list; 

    priority_queue<imp, vector<imp>, cmpimp> q;

    imp s;
    for (int i = 0; i < nImpulses; ++i)
    {
        i_numb = i_numb + 1; // Starts AA index counter at 1, so 0 can be everything in the queue

        s.time = impulse[i];
        s.impind = i_numb;
        s.ind = 0;
        s.pathOg = 1;
        s.pathEnd = 1;
        s.AA_or_not = 1;
        q.push(s);

        s.time = impulse[i];
        s.impind = i_numb;
        s.ind = 10;
        s.pathOg = 2;
        s.pathEnd = 2;
        s.AA_or_not = 1;
        q.push(s);

        
    }
    
    for (int i_q_in = 0; i_q_in < n_q_in/5; ++i_q_in)
    {
        s.time      = q_in[i_q_in*5];
        s.ind       = q_in[i_q_in*5 + 1];
        s.impind    = 0;// Sets AA index to 0 if it comes from the queue, old code: q_in[i_q_in*5 + 2];
        s.pathOg    = q_in[i_q_in*5 + 3];
        s.pathEnd   = q_in[i_q_in*5 + 4];
        s.AA_or_not = 0;
        q.push(s);
    }


    while (!q.empty())
    {        

        // If it get stuck
        numberOfActivations = numberOfActivations + 1;
        if (numberOfActivations > nSave-2)
        {
            plhs[3] = mxCreateCellMatrix(1, 1);
            double* dummy = mxGetPr(plhs[3]);
            dummy = 0;
            break;
        }


        s = q.top();

        if ( (s.AA_or_not == 1)&(s.ind == 0) ) // Saves the current queue if the impulse is coming from the atria
        {
            q_AA_list.push_back(q);

            // To get the States out for each impulse
            for (int i_21 = 0; i_21 < 21; ++i_21) {
                RT_big[(s.impind-1)*21 + i_21] = RT[i_21];
            }
        }

        q.pop();
                
        next_time = s.time;
        next_pos = s.ind;
        pathOg = s.pathOg;
        pathEnd = s.pathEnd;

        s.AA_or_not = 0;

        if (next_time >= RT[next_pos])
        {
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

                if (next_pos==19)
                {
                    s.time = AT[next_pos]+DL[next_pos];
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
                AA_ind_out[c] = s.impind;
                c = c + 1;


                // Switch to new parameters
                for (int i_3 = 0; i_3 < 3; ++i_3) {
                apd_f[i_3] = x_new[i_3];
                apd_s[i_3] = x_new[i_3+3];
                delay_f[i_3] = x_new[i_3+6];
                delay_s[i_3] = x_new[i_3+9];
                }

                
                if (N_RR == c)
                {

                    // Code for outputing the queue for each AA
                    const char *fieldNames[] = {"time"};  // Declare fieldNames here
                    
                    plhs[3] = mxCreateCellMatrix(1, q_AA_list.size());
                    
                    // Convert and store each priority_queue in the cell array
                    for (size_t i = 0; i < q_AA_list.size(); ++i) {
                        const priority_queue<imp, vector<imp>, cmpimp>& q = q_AA_list[i];
                    
                        // Create a matrix for each priority_queue
                        mxArray *qMatrix = mxCreateDoubleMatrix(q.size(), 6, mxREAL);
                        double *qMatrixData = mxGetPr(qMatrix);
                    
                        size_t j = 0;
                    
                        // Create a temporary copy of the priority_queue
                        priority_queue<imp, vector<imp>, cmpimp> tempQ = q;
                    
                        // Use the temporary copy to extract elements without modifying the original
                        while (!tempQ.empty()) {
                            const imp& qElement = tempQ.top();
                    
                            // Populate the matrix with the required information
                            qMatrixData[j] = qElement.time;
                            qMatrixData[j + q.size()] = static_cast<double>(qElement.ind);
                            qMatrixData[j + 2 * q.size()] = static_cast<double>(qElement.impind);
                            qMatrixData[j + 3 * q.size()] = static_cast<double>(qElement.pathOg);
                            qMatrixData[j + 4 * q.size()] = static_cast<double>(qElement.pathEnd);
                            qMatrixData[j + 5 * q.size()] = static_cast<double>(qElement.AA_or_not);
                    
                            tempQ.pop();
                            ++j;

                        }
                    
                        mxSetCell(plhs[3], i, qMatrix);
                    }
                    
                    break;
            }

            }
                        

        
        } // end if
        
        
    } // end while

    if (q.empty())
    {
        plhs[3] = mxCreateCellMatrix(1, 1);
        double* dummy = mxGetPr(plhs[3]);
        dummy = 0;
    }
    
    //mxFree(DL);
    mxFree(AT);

    //mxFree(index);
}


