#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double k = 1.4, km1 = k-1, kp1 = k+1;

vector<double> isentropic(int n, double argv)
{
    double P0_P, T0_T, rho0_rho, M_star, A_A_star, C;
    if(argv<0)
    {
        throw invalid_argument("Invalid input :: input can't be negative\n");
    }
    double M;

    // ======================= 1 calculate mach no. =======================
    if(n==1) // mach no.
    {
        M = argv;
    }
    else if(n==2) // P/Po
    {
        P0_P = argv;
        M = sqrt(2 * (pow(argv, -(km1)/k) - 1) / (km1));
    }
    else if(n==3) // T/To
    {
        T0_T = argv;
        M = sqrt(2 * (1/argv - 1) / (km1));
    }
    else if(n==4) // rho/rho_o
    {
        rho0_rho = argv;
        M = sqrt( 2 * (pow(argv, km1) - 1) / (km1) );
    }
    // the algorithm for area ratio (subsonic & supersonic come from Domingo Tavella, Robust and Efficient Inversion of the Expansion Area Ratio vs. Mach Number Relationship (Stodola’s Equation), November 12, 2014)
    else if(n==5) // A/A* subsonic
    {
        A_A_star = argv;
        if(argv<1)
        {
            throw invalid_argument("Area ratio can't be less than 1\n");
        }
        if(argv>0.9999)
        {
            M = 1.0;
        }
        else
        {
            double f, fd, x, e;
            x = 1+1.0/4*(argv-1)*kp1*(1-sqrt(1+8.0/(argv-1)/kp1));
            // int cnt = 0;
            do {
                M = x;
                f = sgn(1-M)*sqrt(kp1/2.0/km1 * log(2.0/kp1 + km1/kp1*pow(M,2)) - log(M))   -   sqrt(log(argv));
                fd = sgn(1-M)*(M/(2.0/kp1+km1/kp1*pow(M,2)) - 1.0/M)   /   (2*sqrt(kp1/2.0/km1*log(2/kp1+km1/kp1*pow(M,2)) - log(M)));
                x = M - f/fd;
                e = abs(M-x);
                // cnt++;
            } while(e>1e-8);
            // cout << "cnt = " << cnt << "\n";
            M = x;
        }
    }
    else if(n==6) // A/A* supersonic
    {
        A_A_star = argv;
        if(argv<1)
        {
            throw invalid_argument("Area ratio can't be less than 1\n");
        }
        else
        {
            if(argv<1.0001)
            {
                M = 1.0;
            }
            double f, fd, x, e;
            x = pow(argv,km1/2.0)*pow(kp1/km1,kp1/4.0);
            // int cnt = 0;
            do {
                M = x;
                f = sgn(M-1)*sqrt(kp1/2.0/km1 * log(2.0/kp1 + km1/kp1*pow(M,2.0)) - log(M))   -   sqrt(log(argv));
                fd = sgn(M-1)*(M/(2.0/kp1+km1/kp1*pow(M,2.0)) - 1.0/M)   /   (2.0*sqrt(kp1/2.0/km1*log(2.0/kp1+km1/kp1*pow(M,2.0)) - log(M)));
                x = M - f/fd;
                e = abs(M-x);
                // cnt++;
            } while(e>1e-8);
            // cout << "cnt = " << cnt << "\n";
            M = x;
        }
    }
    else if(n==7) // M* = v/a*
    {
        M_star = argv;
        M = sqrt(1 / (kp1/(2+pow(argv,2)) - (km1)/2));
    }
    else if(n==8) // C=v/vmax
    {
        C = argv;
        M = sqrt(2 / (1/pow(argv,2) - 1) / (km1));
    }
    else
    {
        throw invalid_argument("Usage: isentropic(int n, double argv), where n is:\n n=1  argv=M\n n=2  argv=P0/P\n n=3  argv=T/To\n n=4  argv=rho/rho_o\n n=5  argv=A/A*_subsonic\n n=6  argv=A/A*_supersonic\n n=7  argv=M*=v/a*\n n=8  argv=C=v/vmax\n\n");
    }

    // ======================= 2 Calculate all other parameters =======================
    if(n!=2)
    {
        P0_P = pow(1+km1/2.0*pow(M,2), k/km1);
    }
    if(n!=3)
    {
        T0_T = 1+km1/2.0*pow(M,2);
    }
    if(n!=4)
    {
        rho0_rho = pow(1+km1/2.0*pow(M,2), 1.0/km1);
    }
    if(n!=5 || n!=6)
    {
        A_A_star = 1.0/M * pow(2/kp1 + km1/kp1*pow(M,2), kp1/2.0/km1);
    }
    if(n!=7)
    {
        M_star = sqrt(kp1/2.0*pow(M,2)  / (1+km1/2.0*pow(M,2)));
    }
    if(n!=8)
    {
        C = sqrt(km1/2.0*pow(M,2) / (1+km1/2.0*pow(M,2)));
    }
    vector<double> ans = {M, P0_P, T0_T, rho0_rho, A_A_star, M_star, C};
    return ans;
}

vector<double> NSW(double M, double theta = M_PI / 2.0, string direction = "forward")
{
    if(direction=="forward")
    {
        double M1, M2, p2_p1, T2_T1, P02_P01;
        if(M==-1)
        {
            cout << "Usage: NSW(double M1, double theta)\n where theta is in rad.\nOutput: vector<double> ans = {M2, p2_p1, T2_T1, P02_P01}";
        }
        M1 = M * sin(theta);

        M2 = sqrt((2 + km1*pow(M1,2))  /  (2.0*k*pow(M1,2) - km1));
        P02_P01 = pow((2+km1*pow(M2,2)) / (2+km1*pow(M1,2)), k/km1)  *  ((2*k*pow(M1,2)-km1) / kp1);
        p2_p1 = 1 + 2.0*k/kp1*(pow(M1,2)-1);
        T2_T1 = (2+km1*pow(M1,2)) * (2*k*pow(M1,2)-km1) / (pow(kp1*M1,2));
        vector<double> ans = {M2, p2_p1, T2_T1, P02_P01};
        return ans;
    }
    throw invalid_argument("The reversed normal shock problem is not yet implemented in the code\n");
    // implement equation for reverse here
}

double mach_angle(double M, double delta)
{
    // This algorithm gives exact solution and is taken from paper: "T.T. Hartley, R. Brandis, and F. Mossayebi, Exact and Approximate Solutions to the Oblique Shock Equations for Real-Time Applications, NASA Contractor Report 187173, August 1991" 
    // ===== calculate theta =====
    double b = -((pow(M,2)+2)/pow(M,2) + k*pow(sin(delta),2));
    double c = (2*pow(M,2)+1)/pow(M,4) + (pow(kp1,2)/4.0 + km1/pow(M,2))*pow(sin(delta),2);
    double d = -pow(cos(delta),2) / pow(M,4);
    double Q = (3*c-pow(b,2)) / 9.0;
    double R = (9*b*c - 27*d - 2*pow(b,3)) / 54.0;
    double D = pow(Q,3) + pow(R,2);
    if(D>0)
    {
        throw invalid_argument("An error occured: D>0, the deflected shock exists only if D<0, probably this mach no. can't make an attached oblique shockwave with this deflection angle\n");
    }
    double dd = M_PI*(R<0);
    double phi = 1/3.0 * (atan(sqrt(-D) / R) + dd);
    double xw = -b/3.0 - sqrt(-Q) * (cos(phi) - sqrt(3)*sin(phi));
    double xs = -b/3.0 + 2*sqrt(-Q) * cos(phi);
    double theta_strong = atan(sqrt(xs / (1-xs)));
    double theta_weak = atan(sqrt(xw / (1-xw)));
    
    return theta_weak; // choose weak or strong solution
}

vector<double> OSW(double M1, double delta)
{
    double theta = mach_angle(M1, delta);
    // ===== calculate everyting else =====
    vector<double> ans = NSW(M1, theta);
    ans[0] = ans[0] / sin(theta-delta);
    return ans;
}

double prandtl_meyer_angle(int n, double argv)
{
    double M, nu;
    if(n==1) // find nu
    {
        M = argv;
        nu = sqrt(kp1/km1) * atan(sqrt(km1/kp1*(pow(M,2)-1))) - atan(sqrt(pow(M,2)-1));
    }
    else if(n==2) // find M
    {
        // Method 1: ===== Newton-Raphson =====
        // nu = argv;
        // double f, fd, x=1, e;
        // // int cnt = 0;
        // do {
        //     M = x;
        //     f = sqrt(kp1/km1) * atan(sqrt(km1/kp1*(pow(M,2)-1))) - atan(sqrt(pow(M,2)-1)) - nu;
        //     fd = pow(kp1,3/2.0)/(2*km1 + pow(M*km1,2)) - 1.0/pow(M,2);
        //     x = M - f/fd;
        //     e = abs(x-M);
        //     // cnt++;
        // } while(e>1e-9);
        // Method 2: ===== Approximate Inversion of the equation ===== 
        // taken from data sheet Dr. Hesham El-Banna
        double nu_max = M_PI/2.0*(sqrt(6)-1);
        double nu_bar = pow(argv/nu_max, 2/3.0);
        M = (1.0 + 1.3604*nu_bar + 0.0962*pow(nu_bar,2) - 0.5127*pow(nu_bar,3)) / (1.0 - 0.6722*nu_bar - 0.3278*pow(nu_bar,2));
        // cout << "Process took " << cnt << " iterations\n";
    }
    return (n==1)*nu + (n==2)*M;
}

vector<double> EW(double M1, double delta)
{
    double nu1 = prandtl_meyer_angle(1, M1);
    double nu2 = nu1 + delta;
    double M2 = prandtl_meyer_angle(2, nu2);
    vector<double> ans = isentropic(1, M2);
    return ans;
}
