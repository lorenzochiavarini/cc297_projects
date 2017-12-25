#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

// AUXILIARY FUNCTIONS

double distEuclid(double x_a, double y_a, double x_b, double y_b)
/* THIS FUNCTION COMPUTES THE EUCLIDIAN DISTANCE BETWEEN TWO
POINTS IN THE SPACE */
{
    double d;
    d = sqrt((x_b - x_a)*(x_b - x_a) + (y_b - y_a)*(y_b - y_a));
    return d;
}

double x_airfoil(double radius, double angle)
/*THIS FUNCTION COMPUTES THE X COORDINATE OF THE BICONVEX AIRFOIL */
{
    double x_airfoil;
    
    if (angle <= M_PI){
        x_airfoil = radius*cos(angle);
    }
    
    
    return x_airfoil;
}

double y_airfoil(double x_coord, int geometry)
/*THIS FUNCTION COMPUTES THE Y COORDINATE OF THE AIRFOIL */
{
    double y_airfoil;

    if (geometry == 0){
        double t = 0.1;
        y_airfoil = 2.0*t*x_coord*(1.0 - x_coord);
    }

    else {
        double t = 0.12;
        y_airfoil = 5.0*t*(0.2969*pow(x_coord,0.5) - 0.1260*x_coord - 0.3516*pow(x_coord,2.0) \
                    + 0.2843*pow(x_coord,3.0) - 0.1036*pow(x_coord,4.0));
    }
    
    return y_airfoil;
}

int main (void){
    double r0 = 0.5, r1 = 12.0*r0;
    int ksiMax = 93, etaMax = 15;
    int geometry = 1;
    double x_coord[ksiMax][etaMax],y_coord[ksiMax][etaMax];
    double xn[ksiMax][etaMax],xnp1[ksiMax][etaMax];
    double yn[ksiMax][etaMax],ynp1[ksiMax][etaMax];
    double Lx[ksiMax][etaMax],Ly[ksiMax][etaMax];
    double Cx[ksiMax][etaMax],Cy[ksiMax][etaMax];
    double H[ksiMax][etaMax], I[ksiMax][etaMax];
    double a[ksiMax][etaMax], b[ksiMax][etaMax], c[ksiMax][etaMax], d[ksiMax][etaMax];
    double rx=-999999.0,ry=-999999.0;
    double stretching = 1.15;
    double R[etaMax],phiAngle[ksiMax],xRout[ksiMax][1],yRout[ksiMax][1];
    
    FILE *fp;
    
    /* ALBEGRAIC MESH GENERATION */
    
    // definition of concentric circles radii
    for (int j = 0;j < etaMax; j++){
        double p = ((double) j)/((double) etaMax-1);
        R[j] = r0 * (pow((r1/r0),p));
    }
	
    // equally spaced eta lines
    for (int i = 0;i < ksiMax; i++){
	phiAngle[i] = (((double) i)/((double) ksiMax-1)) * 2.0 * M_PI;
	xRout[i][0] = R[etaMax - 1] * cos(phiAngle[i]);
	yRout[i][0] = R[etaMax - 1] * sin(phiAngle[i]);
    }
	
    // x coordinates of each point in the circles
    for (int j = 0;j < etaMax; j++){
        for (int i = 0;i < ksiMax; i++){
	        x_coord[i][j] = R[j] * cos(phiAngle[i]);
	    }
    }
    
    // y coordinates of each point + airfoil wall
    for (int j = 0;j < etaMax; j++){
        for (int i = 0;i < ksiMax; i++){
            double x_airfoil = fabs(x_coord[i][j] + r0);
            if (j == 0){
                if (phiAngle[i] < M_PI){
                    y_coord[i][j] = y_airfoil(x_airfoil,geometry);
                }
                else {
                   y_coord[i][j] = -1.0*y_airfoil(x_airfoil,geometry);
                }
            }
            else {
                y_coord[i][j] = R[j] * sin(phiAngle[i]);
            }
        }
    }
    
    // finding x and y coordinates using Euclidian distance
    for (int i = 0;i < ksiMax; i++){
	    double b = 0.0;
	    for (int j = 1;j < etaMax-1; j++){
            double d = distEuclid(x_coord[i][0],y_coord[i][0],xRout[i][0],yRout[i][0]);
            double a = d*(stretching-1.0)/(pow(stretching,((double) etaMax-1)));
            b = (b + a)*stretching;
            x_coord[i][j] = (b/d)*(xRout[i][0]-x_coord[i][0])+x_coord[i][0];
            y_coord[i][j] = (b/d)*(yRout[i][0]-y_coord[i][0])+y_coord[i][0];
        }
    }

    // WRITTING ALGEBRAIC MESH FILE

    fp = fopen("algebraicMesh.dat","w+");
    
    fprintf(fp,"TITLE = \"Projeto 02\"\nVARIABLES = \"X\",\"Y\",\"PHI\",\"R\"\nZONE T = \"Zone-one\",I=%d, J=%d, F=POINT\n",ksiMax,etaMax);
    for (int j = 0;j < etaMax; j++){
        for (int i = 0;i < ksiMax; i++){
	        fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",x_coord[i][j],y_coord[i][j],0.0,0.0);
	    }
	}

    fclose(fp);
    // ----------------------------------------------------------------
    // ----------------------------------------------------------------
    // ----------------------------------------------------------------
    // ----------------------------------------------------------------

    // initializing matrices and vectors for iterations
    for (int j = 0;j < etaMax; j++){
        for (int i = 0;i < ksiMax; i++){
	        xn[i][j] = x_coord[i][j];
	        yn[i][j] = y_coord[i][j];
	        Cx[i][j] = 0.0;
	        Cy[i][j] = 0.0;
                Lx[i][j] = 0.0;
                Ly[i][j] = 0.0;
                H[i][j] = 0.0;
                I[i][j] = 0.0;
                a[i][j] = 0.0;
                b[i][j] = 0.0;
                c[i][j] = 0.0;
        }
    }
    
    // NUMERICAL SCHEME
    int it = 0,itMax = 5;
    double tol = 1e-20;
    double alpha = 1.0, omega = 1.2;

    fp = fopen("residuals.dat","w+");
    do{
        // iterating for x

        for (int j = 1;j < etaMax-1; j++){
            for(int i = 1; i < ksiMax-1; i++){
                a[i][j] = 0.25*(xn[i][j+1]-xn[i][j-1])*(xn[i][j+1]-xn[i][j-1])+0.25*(yn[i][j+1]-yn[i][j-1])*(yn[i][j+1]-yn[i][j-1]);
                b[i][j] = 0.25*(xn[i+1][j]-xn[i-1][j])*(xn[i][j+1]-xn[i][j-1])+0.25*(yn[i+1][j]-yn[i-1][j])*(yn[i][j+1]-yn[i][j-1]);
                c[i][j] = 0.25*(xn[i+1][j]-xn[i-1][j])*(xn[i+1][j]-xn[i-1][j])+0.25*(yn[i+1][j]-yn[i-1][j])*(yn[i+1][j]-yn[i-1][j]);
            }
        }

        for (int j = 1;j < etaMax-1; j++){
            for(int i = 1; i < ksiMax-1; i++){
                double xee = (xn[i+1][j] - 2.0*xn[i][j] + xn[i-1][j]);
                double xen = 0.25*(xn[i+1][j+1] - xn[i+1][j-1] - xn[i-1][j+1] + xn[i-1][j-1]);
                double xnn = (xn[i][j+1] - 2.0*xn[i][j] + xn[i][j-1]);
                Lx[i][j] = a[i][j]*xee -2.0*b[i][j]*xen + c[i][j]*xnn ;
                if (Lx[i][j] > rx){
                    rx = Lx[i][j];
                }
            }
        }

        // first step for ADI
        for (int j = 1;j < etaMax-1; j++){
            for(int i = 1; i < ksiMax-1; i++){
                double N = -(1.0/alpha)*(alpha + 2.0*a[i][j]);
                H[i][j] = -(omega*Lx[i][j] + (a[i][j]/alpha)*H[i+1][j] + (a[i][j]/alpha)*H[i-1][j])/N;
            }
        }

        // second step for ADI
        for (int j = 1;j < etaMax-1; j++){
            for(int i = 1; i < ksiMax-1; i++){
                double N = (alpha - 2.0*c[i][j]);
                Cx[i][j] = (H[i][j] + c[i][j]*Cx[i][j+1] + c[i][j]*Cx[i][j-1])/N;
            }
        }

        for (int j = 1;j < etaMax-1; j++){
            for(int i = 1; i < ksiMax-1; i++){
                xnp1[i][j] = xn[i][j] + Cx[i][j];
            }
        }

        for (int j = 1;j < etaMax-1; j++){
            for(int i = 1; i < ksiMax-1; i++){
                xn[i][j] = xnp1[i][j];
            }
        }
        

        // iterating for y

        for (int j = 1;j < etaMax-1; j++){
            for(int i = 1; i < ksiMax-1; i++){
                double yee = (yn[i+1][j] - 2.0*yn[i][j] + yn[i-1][j]);
                double yen = 0.25*(yn[i+1][j+1] - yn[i+1][j-1] - yn[i-1][j+1] + yn[i-1][j-1]);
                double ynn = (yn[i][j+1] - 2.0*yn[i][j] + yn[i][j-1]);
                Ly[i][j] = a[i][j]*yee -2.0*b[i][j]*yen + c[i][j]*ynn ;               
                if (Ly[i][j] > ry){
                    ry = Ly[i][j];
                }
            }
        }

        // first step for ADI
        for (int j = 1;j < etaMax-1; j++){
            for(int i = 1; i < ksiMax-1; i++){
                double N = -(1.0/alpha)*(alpha + 2.0*a[i][j]);
                I[i][j] = -(omega*Ly[i][j] + (a[i][j]/alpha)*I[i+1][j] + (a[i][j]/alpha)*I[i-1][j])/N;
            }
        }

        // second step for ADI
        for (int j = 1;j < etaMax-1; j++){
            for(int i = 1; i < ksiMax-1; i++){
                double N = (alpha - 2.0*c[i][j]);
                Cy[i][j] = (I[i][j] + c[i][j]*Cy[i][j+1] + c[i][j]*Cy[i][j-1])/N;
            }
        }

        for (int j = 1;j < etaMax-1; j++){
            for(int i = 1; i < ksiMax-1; i++){
                ynp1[i][j] = yn[i][j] + Cy[i][j];
            }
        }

        for (int j = 1;j < etaMax-1; j++){
            for(int i = 1; i < ksiMax-1; i++){
                yn[i][j] = ynp1[i][j];
            }
        }

        fprintf(fp,"%d\t%lf\n",it,log10(rx));

        it++;
    
    }while (it <= itMax && ry > tol);
    fclose(fp);

    // WRITTING ITERATED MESH DATA FILE
    fp = fopen("iteratedmsh.dat","w+");
    fprintf(fp,"TITLE = \"Projeto 02\"\nVARIABLES = \"X\",\"Y\",\"PHI\",\"R\"\nZONE T = \"Zone-one\",I=%d, J=%d, F=POINT\n",ksiMax,etaMax);
    for (int j = 0;j < etaMax; j++){
	    for (int i = 0;i < ksiMax; i++){
	        fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",x_coord[i][j],y_coord[i][j],0.0,0.0);
	    }
	}
    fclose(fp);
    
    return 0;
}
