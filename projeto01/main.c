#include<stdlib.h>
#include<stdio.h>
#include<math.h>

// AUXILIARY FUNCTIONS
double boundary_phi(double u,double x)
/* FUNCTION TO CALCULATE THE BOUNDARY
CONDITIONS AT THE TOP AND AT BOTH RIGHT
AND LEFT SIDES OF THE DOMAIN*/
{
    double phi;
    phi = u * x;
    return phi;
}

double airfoil_phi(double u,double t,double x)
/* FUNCTION TO CALCULATE THE BOUNDARY
CONDITIONS AT THE BOTTOM OF THE DOMAIN,
AT THE AIRFOIL LOCATION */
{
    double phi;
    phi = u * (2.0*t - 4.0*t*x);
    return phi;
}

double calcVelocityU(double phi_ip1j, double phi_ij, double phi_im1j, double x_ip1j, double x_ij, double x_im1j, int i, int IMAX)
/* THIS FUNCTION CALCULATES THE X COMPONENT OF VELOCITY VECTOR USING FINITE DIFFERENCE SCHEME */
{
    double u,dx;
    if (i == 0) // FORWARD DIFFERENCE SCHEME FOR THE LEFT BOUNDARY
    {
        dx = x_ip1j-x_ij;
        u = (phi_ip1j - phi_ij)/dx;
    }
    else if ( i > 0 && i < IMAX-1) // CENTERED DIFFERENCE SCHEME FOR INTERNAL NODES
    {
        dx = x_ip1j-x_im1j;
        u = (phi_ip1j - phi_im1j)/dx;
    }
    else if (i == IMAX-1) // BACKWARD DIFFERENCE SCHEME FOR THE LEFT BOUNDARY
    {
        dx = x_ij-x_im1j;
        u = (phi_ij - phi_im1j)/dx;
    }
    return u;
}

double calcVelocityV(double phi_ijp1, double phi_ij, double phi_ijm1, double y_ijp1, double y_ij, double y_ijm1, int j, int JMAX)
/* THIS FUNCTION CALCULATES THE Y COMPONENT OF VELOCITY VECTOR USING FINITE DIFFERENCE SCHEME */
{
    double v,dy;

    if (j == 0) // FORWARD DIFFERENCE SCHEME FOR THE LOWER BOUNDARY
    {
        dy = y_ijp1-y_ij;
        v = (phi_ijp1 - phi_ij)/dy;
    }
    else if ( j > 0 && j < JMAX-1) // CENTERED DIFFERENCE SCHEME FOR THE INTERNAL NODES
    {
        dy = y_ijp1-y_ijm1;
        v = (phi_ijp1 - phi_ijm1)/dy;
    }
    else if (j == JMAX-1) // BACKWARD DIFFERENCE SCHEME FOR THE UPPER BOUNDARY
    {
        dy = y_ij-y_ijm1;
        v = (phi_ij - phi_ijm1)/dy;
    }
    return v;
}

double calcCpIncomp(double u, double u_inf)
/* THIS FUNCTION CALCULATES THE PRESSURE COEFFICIENT (Cp) FOR AN INCOMPRESSIBLE FLOW */
{
    double cp;
    cp = 1.0 - ((u*u)/(u_inf*u_inf));
    return cp;
}

double modU(double u, double v, double w)
/* THIS FUNCTION CALCULATES THE MAGNITUDE OF VECTOR U */
{
    double U;
    U = sqrt(u*u + v*v + w*w);
    return U;
}

int main( )
/* THIS IS THE MAIN CODE. HERE IS WHERE ONE SETS THE PARAMETERS FOR THE SIMULATIONS
AS WELL AS SELECT THE CASE TO RUN */
{
    int write_output = 1, write_cp = 0, write_mesh = 0;     // write_file = 1 - yes, 2 - no
    int i,j,n=0;
    int ILE = 11, ITE = 31, IMAX = 41, JMAX = 12;
    int caso = 1;                                           // 1 -> t = 0.05, u_inf = 1.0; 2 -> t = 0.10, u_inf = 1.0
    int method = 2;                                         // 1 -> Point-Jacobi, 2 -> Gauss-Seidel, 3 -> SOR, 4 -> Line Gauss-Seidel, 5 -> SLOR

    int k = 15000;                                          // iteration parameters -> max iterations (k), tolerance (tol)
    double tol = 1.0e-20;

    double dx, dy, XSF = 1.25, YSF = 1.25;                  // dx, dy -> mesh spacing in x and y directions, respectively
    double u_inf, t;                                        // u_inf and t are input data
    double s = 0.0,err=0.0,res,beta;                        // s, err, res and beta are support variables
    double x[IMAX][JMAX],y[IMAX][JMAX];                     // arrays with x and y coordinates for each i,j point in the grid
    double U[IMAX][JMAX], u[IMAX][JMAX],v[IMAX][JMAX];      // u and v components of velocity
    double cp[IMAX][JMAX];                                  // array with cp distribution
    double Lphi[IMAX][JMAX], Cphi[IMAX][JMAX];              // Lphi is the residual operator; cp is the pressure coefficient; Cphi is the difference between the solution at n+1 and n
    double phi[IMAX][JMAX],phi_np1[IMAX][JMAX];

    FILE *fp;

    // selecting the case to run
    if (caso == 1)
    {
        u_inf = 1.0;
        t = 0.05;
    }
    else if(caso == 2)
    {
        u_inf = 1.0;
        t = 0.10;
    }
    // initializing variables
    for (i = 0; i < IMAX; i++)
    {
        for (j = 0; j < JMAX; j++)
        {
            x[i][j] = -999999;
            y[i][j] = -999999;
            Lphi[i][j] = -999999;
            Cphi[i][j] = -999999;
            phi[i][j] = -999999;
            phi_np1[i][j] = -999999;
            cp[i][j] = -999999;
            u[i][j] = -999999;
            v[i][j] = -999999;
            U[i][j] = -999999;
        }
    }
    // ===============================================================
    // generating computational grid
    // ===============================================================
    dx = 1.0 / (double)(ITE - ILE); // as defined in the proposed work
    // --------------------- y coordinates -------------------------
    for (i = 0; i < IMAX; i++)
    {
        for (j = 0; j < JMAX; j++)
        {
            if (j == 0)
                y[i][j] = -dx/2.0;
            else if (j == 1)
                y[i][j] = dx/2.0;
            else
                y[i][j] = y[i][j-1] + (y[i][j-1] - y[i][j-2])*YSF;
        }
    }
    // --------------------- x coordinates -------------------------
    // x coordinates for nodes between leading and trailing edges
    for (i = ILE - 1; i < ITE; i++)
    {
        for (j = 0; j < JMAX; j++)
        {
            x[i][j] = (double)(i - (ILE-1))*dx;
        }
    }
    // x coordinates for nodes after the trailing edge
    for (i = ITE; i < IMAX; i++)
    {
        for (j = 0; j < JMAX; j++)
        {
            x[i][j] = x[i-1][j] + (x[i-1][j] - x[i-2][j])*XSF;
        }
    }
    // x coordinates for nodes prior to the leading edge
    for (i = ILE - 2; i != -1; i-=1)
    {
        for (j = 0; j < JMAX; j++)
        {
            x[i][j] = x[i+1][j] + (x[i+1][j] - x[i+2][j])*XSF;
        }
    }
    // --------------------------------------------------------------
    // --------------------------------------------------------------
    // INITIALIZING PHI
    for (i = 0;i < IMAX; i++)
    {
        for (j = 0; j < JMAX; j++)
        {
            phi[i][j] = boundary_phi(u_inf,x[i][j]);
        }
    }
    for (i = 0; i < IMAX; i++)
    {
        if (x[i][0] >= 0.0 && x[i][0] <= 1.0)
        {
            phi[i][0] = phi[i][1] - (y[i][1] - y[i][0])*airfoil_phi(u_inf,t,x[i][0]);
        }
        else
            phi[i][0] = phi[i][1];
    }

    // --------------------------------------------------------------
    // --------------------------------------------------------------
    // starting phi_np1 equal to phi arrays
    for (i = 0; i < IMAX; i++)
    {
        for (j = 0; j < JMAX; j++)
        {
            phi_np1[i][j] = phi[i][j];
            Cphi[i][j] = 0;
        }
    }
    // --------------------------------------------------------------
    // --------------------------------------------------------------
    // POINT JACOBI
    if (method == 1)
    {
        printf("\n___________________________________________________\n");
        printf("         >>> RUNNING POINT JACOBI METHOD <<<         \n");
        printf("\n___________________________________________________\n");
        fp = fopen("PJresidual.dat","w+");
        do
        {
            res = 0.0;
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)
                {
//                    dx = (x[i+1][j] - x[i-1][j])/2.0;
//                    dy = (y[i][j+1] - y[i][j-1])/2.0;
                    // CALCULATING THE RESIDUAL
                    Lphi[i][j] = (2.0/(x[i+1][j]-x[i-1][j]))*(((phi[i+1][j]-phi[i][j])/(x[i+1][j]-x[i][j]))-((phi[i][j]-phi[i-1][j])/(x[i][j]-x[i-1][j])))
                                +(2.0/(y[i][j+1]-y[i][j-1]))*(((phi[i][j+1]-phi[i][j])/(y[i][j+1]-y[i][j]))-((phi[i][j]-phi[i][j-1])/(y[i][j]-y[i][j-1])));
                    // GETTING THE MAXIMUM RESIDUAL
                    if (fabs(Lphi[i][j]) > res)
                    {
                        res = fabs(Lphi[i][j]);
                    }
                }
            }
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)
                {
                    dx = (x[i+1][j] - x[i-1][j])/2.0;
                    dy = (y[i][j+1] - y[i][j-1])/2.0;
                    // SOLVING FOR Cphi
                    Cphi[i][j] = -(1.0/((-2.0/pow(dx,2.0))+(-2.0/pow(dy,2.0))))*Lphi[i][j];
                }
            }
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)
                {
                    // CALCULATING THE POTENTIAL AT n+1
                    phi_np1[i][j] = phi[i][j] + Cphi[i][j];
                }
            }
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)
                {
                    // UPDATING phi
                    phi[i][j] = phi_np1[i][j];
                }
            }
            // UPDATING FOR j = 0
            for (i = 0; i < IMAX; i++)
            {
                if (x[i][0] >= 0.0 && x[i][0] <= 1.0)
                {
                    phi[i][0] = phi[i][1] - (y[i][1] - y[i][0])*airfoil_phi(u_inf,t,x[i][0]);
                }
                else
                    phi[i][0] = phi[i][1];
            }
            fprintf(fp,"%d\t %f\n",n,log10(res));
            n++;
        }while (n<=k && res > tol);
        fclose(fp);
    }

    // POINT GAUSS-SEIDEL METHOD
    else if (method == 2)
    {
        printf("\n___________________________________________________\n");
        printf("      >>> RUNNING POINT GAUSS-SEIDEL METHOD <<<      \n");
        printf("\n___________________________________________________\n");
        fp = fopen("PGSresidual.dat","w+");
        do
        {
            res = 0.0;
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)
                {
                    // CALCULATING THE RESIDUAL
                    Lphi[i][j] = (2.0/(x[i+1][j]-x[i-1][j]))*(((phi[i+1][j]-phi[i][j])/(x[i+1][j]-x[i][j]))-((phi[i][j]-phi[i-1][j])/(x[i][j]-x[i-1][j])))+\
                                +(2.0/(y[i][j+1]-y[i][j-1]))*(((phi[i][j+1]-phi[i][j])/(y[i][j+1]-y[i][j]))-((phi[i][j]-phi[i][j-1])/(y[i][j]-y[i][j-1])));
                    // GETTING THE MAXIMUM RESIDUAL
                    if (fabs(Lphi[i][j]) > res)
                    {
                        res = fabs(Lphi[i][j]);
                    }
                }
            }
                    
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)

                {
                    dx = (x[i+1][j] - x[i-1][j])/2;
                    dy = (y[i][j+1] - y[i][j-1])/2;
                    // CALCULATING FOR Cphi
                    Cphi[i][j] = -(1.0/((-2.0/pow(dx,2.0))+(-2.0/pow(dy,2.0))))*((Cphi[i-1][j]/pow(dx,2.0)) + (Cphi[i][j-1]/pow(dy,2.0)) + Lphi[i][j]);
                }
            }
                 
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)

                {
                    // CALCULATING POTENTIAL AT n+1
                    phi_np1[i][j] = phi[i][j] + Cphi[i][j];
                }
            }
     
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)

                {
                    // UPDATING PHI
                    phi[i][j] = phi_np1[i][j];
                }
            }
     
            // UPDATING PHI FOR j=0
            for (i = 0; i < IMAX; i++)
            {
                if (x[i][0] >= 0.0 && x[i][0] <= 1.0)
                {
                    phi[i][0] = phi[i][1] - (y[i][1] - y[i][0])*airfoil_phi(u_inf,t,x[i][0]);
                }
                else
                    phi[i][0] = phi[i][1];
            }
            fprintf(fp,"%d\t %f\n",n,log10(res));
            n++;
        }while (n<=k && res > tol);
        fclose(fp);
    }
    
    // SUCCESSIVE OVERRELAXATION METHOD
    else if (method == 3)
    {
        printf("\n___________________________________________________\n");
        printf("   >> RUNNING SUCCESSIVE OVERRELAXATION METHOD <<<   \n");
        printf("\n___________________________________________________\n");
        fp = fopen("SORresidual.dat","w+");
        double r = 1.85;
        do
        {
            res = 0.0;
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)
                {
                    // CALCULATING THE RESIDUAL
                    Lphi[i][j] = (2.0/(x[i+1][j]-x[i-1][j]))*(((phi[i+1][j]-phi[i][j])/(x[i+1][j]-x[i][j]))-((phi[i][j]-phi[i-1][j])/(x[i][j]-x[i-1][j])))+\
                                +(2.0/(y[i][j+1]-y[i][j-1]))*(((phi[i][j+1]-phi[i][j])/(y[i][j+1]-y[i][j]))-((phi[i][j]-phi[i][j-1])/(y[i][j]-y[i][j-1])));
                    // GETTING THE MAXIMUM RESIDUAL
                    if (fabs(Lphi[i][j]) > res)
                    {

                        res = fabs(Lphi[i][j]);
                    }
                }
            }
                    
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)

                {
                    dx = (x[i+1][j] - x[i-1][j])/2;
                    dy = (y[i][j+1] - y[i][j-1])/2;
                    // CALCULATING FOR Cphi
                    Cphi[i][j] = -(r/((-2.0/pow(dx,2.0))+(-2.0/pow(dy,2.0))))*((Cphi[i-1][j]/pow(dx,2.0)) + (Cphi[i][j-1]/pow(dy,2.0)) + Lphi[i][j]);
                }
            }
                 
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)

                {
                    // CALCULATING POTENTIAL AT n+1
                    phi_np1[i][j] = phi[i][j] + Cphi[i][j];
                }
            }
     
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)

                {
                    // UPDATING PHI
                    phi[i][j] = phi_np1[i][j];
                }
            }
     
            // UPDATING PHI FOR j=0
            for (i = 0; i < IMAX; i++)
            {
                if (x[i][0] >= 0.0 && x[i][0] <= 1.0)
                {
                    phi[i][0] = phi[i][1] - (y[i][1] - y[i][0])*airfoil_phi(u_inf,t,x[i][0]);
                }
                else
                    phi[i][0] = phi[i][1];
            }
            fprintf(fp,"%d\t %f\n",n,log10(res));
            n++;
        }while (n<=k && res > tol);
        fclose(fp);
    }

    // LINE GAUSS-SEIDEL METHOD
    else if (method == 4)
    {
        printf("\n___________________________________________________\n");
        printf("       >> RUNNING LINE GAUSS-SEIDEL METHOD <<<       \n");
        printf("\n___________________________________________________\n");
        fp = fopen("LGSresidual.dat","w+");
        do
        {
            res = 0.0;
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)
                {
                    // CALCULATING THE RESIDUAL
                    Lphi[i][j] = (2.0/(x[i+1][j]-x[i-1][j]))*(((phi[i+1][j]-phi[i][j])/(x[i+1][j]-x[i][j]))-((phi[i][j]-phi[i-1][j])/(x[i][j]-x[i-1][j])))+\
                                +(2.0/(y[i][j+1]-y[i][j-1]))*(((phi[i][j+1]-phi[i][j])/(y[i][j+1]-y[i][j]))-((phi[i][j]-phi[i][j-1])/(y[i][j]-y[i][j-1])));
                    // GETTING THE MAXIMUM RESIDUAL
                    if (fabs(Lphi[i][j]) > res)
                    {

                        res = fabs(Lphi[i][j]);
                    }
                }
            }
                    
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)

                {
                    dx = (x[i+1][j] - x[i-1][j])/2;
                    dy = (y[i][j+1] - y[i][j-1])/2;
                    // CALCULATING FOR Cphi
                    Cphi[i][j] = -(1.0/((-2.0/pow(dx,2.0))+(-1.0/(dy*(y[i][j+1]-y[i][j])))+(-1.0/(dy*(y[i][j]-y[i][j-1])))))*((Cphi[i-1][j]/pow(dx,2.0)) 
                                 +(Cphi[i][j+1]/(dy*(y[i][j+1]-y[i][j]))) + (Cphi[i][j-1]/(dy*(y[i][j]-y[i][j-1]))) + Lphi[i][j]);
                }
            }
                 
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)

                {
                    // CALCULATING POTENTIAL AT n+1
                    phi_np1[i][j] = phi[i][j] + Cphi[i][j];
                }
            }
     
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)

                {
                    // UPDATING PHI
                    phi[i][j] = phi_np1[i][j];
                }
            }
     
            // UPDATING PHI FOR j=0
            for (i = 0; i < IMAX; i++)
            {
                if (x[i][0] >= 0.0 && x[i][0] <= 1.0)
                {
                    phi[i][0] = phi[i][1] - (y[i][1] - y[i][0])*airfoil_phi(u_inf,t,x[i][0]);
                }
                else
                    phi[i][0] = phi[i][1];
            }
            fprintf(fp,"%d\t %f\n",n,log10(res));

            n++;
        }while (n<=k && res > tol);
        fclose(fp);
    }

    // SUCCESSIVE LINE OVERRELAXATION METHOD
    else if (method == 5)
    {
        printf("\n___________________________________________________\n");
        printf(">>> RUNNING SUCCESSIVE LINE OVERRELAXATION METHOD <<<\n");
        printf("\n___________________________________________________\n");
        fp = fopen("SLORresidual.dat","w+");
        double r = 1.85;
        do
        {
            res = 0.0;
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)
                {
                    // CALCULATING THE RESIDUAL
                    Lphi[i][j] = (2.0/(x[i+1][j]-x[i-1][j]))*(((phi[i+1][j]-phi[i][j])/(x[i+1][j]-x[i][j]))-((phi[i][j]-phi[i-1][j])/(x[i][j]-x[i-1][j])))+\
                                +(2.0/(y[i][j+1]-y[i][j-1]))*(((phi[i][j+1]-phi[i][j])/(y[i][j+1]-y[i][j]))-((phi[i][j]-phi[i][j-1])/(y[i][j]-y[i][j-1])));
                    // GETTING THE MAXIMUM RESIDUAL
                    if (fabs(Lphi[i][j]) > res)
                    {

                        res = fabs(Lphi[i][j]);
                    }
                }
            }
                    
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)

                {
                    dx = (x[i+1][j] - x[i-1][j])/2;
                    dy = (y[i][j+1] - y[i][j-1])/2;
                    // CALCULATING FOR Cphi
//                    Cphi[i][j] = -(1.0/((-2.0/pow(dx,2.0))+(-2.0/pow(dy,2.0))))*(r*((Cphi[i-1][j]/pow(dx,2.0)) + Lphi[i][j]) + (Cphi[i][j-1]/pow(dy,2.0))
//                                 +(Cphi[i][j+1]/pow(dy,2.0)));
                    Cphi[i][j] = -(1.0/((-2.0/pow(dx,2.0))+(-1.0/(dy*(y[i][j+1]-y[i][j])))+(-1.0/(dy*(y[i][j]-y[i][j-1])))))*(r*((Cphi[i-1][j]/pow(dx,2.0)) +Lphi[i][j]) + (Cphi[i][j+1]/(dy*(y[i][j+1]-y[i][j]))) + (Cphi[i][j-1]/(dy*(y[i][j]-y[i][j-1]))));
                }
            }
                 
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)

                {
                    // CALCULATING POTENTIAL AT n+1
                    phi_np1[i][j] = phi[i][j] + Cphi[i][j];
                }
            }
     
            for (i = 1; i < IMAX-1; i++)
            {
                for (j = 1; j < JMAX-1; j++)

                {
                    // UPDATING PHI
                    phi[i][j] = phi_np1[i][j];
                }
            }
     
            // UPDATING PHI FOR j=0
            for (i = 0; i < IMAX; i++)
            {
                if (x[i][0] >= 0.0 && x[i][0] <= 1.0)
                {
                    phi[i][0] = phi[i][1] - (y[i][1] - y[i][0])*airfoil_phi(u_inf,t,x[i][0]);
                }
                else
                    phi[i][0] = phi[i][1];
            }
            fprintf(fp,"%d\t %f\n",n,log10(res));

            n++;
        }while (n<=k && res > tol);
        fclose(fp);
    }

    // CALCULATING VELOCITY COMPONENTS AND PRESSURE COEFFICIENT
    for (i = 0; i < IMAX; i++)
    {
        for (j = 0; j < JMAX; j++)
        {
            u[i][j] = calcVelocityU(phi[i+1][j],phi[i][j],phi[i-1][j],x[i+1][j],x[i][j],x[i-1][j],i,IMAX);
            v[i][j] = calcVelocityV(phi[i][j+1],phi[i][j],phi[i][j-1],y[i][j+1],y[i][j],y[i][j-1],j,JMAX);
        }
    }
    for (i = 0; i < IMAX; i++)
    {
        for (j = 0; j < JMAX; j++)
        {
            U[i][j] = modU(u[i][j],v[i][j],0.0);
            cp[i][j] = calcCpIncomp(U[i][j],u_inf);
        }
    }

    // WRITING CP.DAT FILE
    if (write_cp == 1)
    {
        fp = fopen("cp.dat","w+");
        for (i = 0; i < IMAX; i++)
        {
            if (x[i][0] >= 0.0 && x[i][0] <= 1.0)
            {
                double cpm = 0.5*(cp[i][0]+cp[i][1]);
                fprintf(fp,"%f\t %f\n",x[i][0],-cpm);
            }
        }
        fclose(fp);
    }

    // WRITING OUTPUT.DAT FILE
    if (write_output == 1)
    {
        fp = fopen("output.dat","w+");
        fprintf(fp,"TITLE = \"Projeto 01\"\nVARIABLES = \"X\",\"Y\",\"PHI\",\"u_speed\",\"v_speed\",\"cp\"\nZONE T =\"Zone-one\", I= %d, J= %d, F=POINT\n",IMAX,JMAX);

        for (j = 0; j < JMAX; j++)
        {
            for (i = 0; i < IMAX; i++)
            {
                fprintf(fp,"%f\t %f\t %f\t %f\t %f\t %f\n",x[i][j],y[i][j],phi[i][j],u[i][j],v[i][j],-cp[i][j]);
            }
        }
        fclose(fp);
    }

    // WRITING MSH.DAT FILE
    if (write_mesh == 1)
    {
        fp = fopen("msh.dat","w+");
        fprintf(fp,"TITLE = \"Projeto 01\"\nVARIABLES = \"X\",\"Y\",\"PHI\",\"R\"\nZONE T =\"Zone-one\", I= %d, J= %d, F=POINT\n",IMAX,JMAX);

        for (j = 0; j < JMAX; j++)
        {
            for (i = 0; i < IMAX; i++)
            {
                fprintf(fp,"%f\t %f\t %f\t %f\n",x[i][j],y[i][j],phi_np1[i][j],Lphi[i][j]);
//                fprintf(fp,"%f\t %f\n",x[i][j],y[i][j]);
            }
        }
        fclose(fp);
    }
    return 1;
}
