// This is code to solve a supersonic airfoil problem using expansion theory

#include "Gas_Dynamics.cpp"
using namespace std;

void print_state(int n, vector<double> vec)
{
    vector<string> str;
    if (n == 1) // isentropic
    {
        str = {"M", "P0_P", "T0_T", "rho0_rho", "A_A*", "M*", "C"};
    }
    else if (n == 2) // shock
    {
        str = {"M", "P02_P01", "p2_p1*", "T02_T01"};
    }
    int sz = (int)vec.size();

    cout << "[";
    for (int i = 0; i < sz; i++)
    {
        cout << str[i] << " = " << vec[i];
        if (i != sz - 1)
        {
            cout << ", ";
        }
    }
    cout << "]\n";
}

int main()
{
    // define the problem
    double M_inf = 2, alpha = 5*M_PI/180.0, epsilon = 10*M_PI/180.0; // epsilon is the half wedge angle
    
    cout << "inf --->1 : OSW\n";
    vector<double> state1 = OSW(M_inf, epsilon-alpha);
    double M1 = state1[0];
    double p1_p_inf = state1[1];
    cout << " M1 = " << M1 << ", p1_p_inf = " << p1_p_inf << "\n\n";
    
    cout << "1 ---> 2 : EW\n";
    vector<double> state2 = EW(M1, 2*epsilon);
    double M2 = state2[0];
    double p2_p1 = 1.0/state2[1] * isentropic(1, M1)[1]; // p2_p1 = p2_p0 / p0_p1
    double p2_p_inf = p2_p1 * p1_p_inf;
    cout << " M2 = " << M2 << ", p2_p_inf = " << p2_p_inf << "\n\n";


    cout << "inf --->3 : OSW\n";
    vector<double> state3 = OSW(M_inf, epsilon+alpha);
    double M3 = state3[0];
    double p3_p_inf = state3[1];
    cout << " M3 = " << M3 << ", p3_p_inf = " << p3_p_inf << "\n\n";
    
    cout << "1 ---> 4 : EW\n";
    vector<double> state4 = EW(M3, 2*epsilon);
    double M4 = state4[0];
    double p4_p3 = 1.0/state4[1] * isentropic(1, M3)[1];
    double p4_p_inf = p4_p3 * p3_p_inf;
    cout << " M4 = " << M4 << ", p4_p_inf = " << p4_p_inf << "\n\n";

    // calculate cn, ca, cl, cd
    double cn = 1.0/k/pow(M_inf, 2)*(-p1_p_inf-p2_p_inf+p3_p_inf+p4_p_inf);
    double ca = 1.0/k/pow(M_inf, 2)*(p1_p_inf-p2_p_inf+p3_p_inf-p4_p_inf)*tan(epsilon);
    double cl = cn*cos(alpha) - ca*sin(alpha);
    double cd = ca*cos(alpha) + cn*sin(alpha);
    double cm_le = 1.0/4.0/k/pow(M_inf, 2) * ((p1_p_inf+3*p2_p_inf-p3_p_inf-3*p4_p_inf) + pow(tan(epsilon),2)*(p1_p_inf-p2_p_inf-p3_p_inf+p4_p_inf));
    double cm_c2 = cm_le + 1/2.0*cl;
    cout << "cn = " << cn << ", ca = " << ca << ", cl = " << cl << ", cd = " << cd << ", cm_le = " << cm_le << ", cm_c2 = " << cm_c2 << "\n";
    return 0;
}
