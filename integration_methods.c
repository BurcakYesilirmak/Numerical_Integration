    #include<stdio.h>
    #include<math.h>

    double norm(double r[3]) {   
        double norm_r = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
        return (norm_r);
    }

    double energies(double r[3],double V[3]){
        double norm(double r[]);
        double norm_r=norm(r);
        double v2 =V[0]*V[0]+V[1]*V[1]+V[2]*V[2];
        double ek=0.5*v2;
        double epot =-1/norm_r;
        double conservation= ek + epot;
        return(conservation);
    } 
     /* void accelaration(double r[3], double *a){
    double r2 =r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
    double r3=r2*sqrt(r2);
    for(int i=0;i<3;i++)
    a[i] = -r[i]/r3;
    return;
    } */
    void euler(double r[3],double V[3], double a[3],double vtemp[3]){
         double eo =energies(r,V);
         double dt;
         double norm(double r[]);
         double norm_r = norm(r);
         double r_denom = pow(norm_r,3); 
        a[0] = -r[0]/r_denom;  // x ve a yı değiştir
        a[1] = -r[1]/r_denom;
        a[2] = -r[2]/r_denom;

        r[0] += V[0]*dt;
        r[1] += V[1]*dt;
        r[2] += V[2]*dt;

        V[0] += a[0]*dt;   // eğer V ve r yi yerdeğiştirirsen korunum olan bir kodda
        V[1] += a[1]*dt;   // kapalı bir eğri elde edersin.
        V[2] += a[2]*dt;
        return;
    }

    void leapfrog(double *r[3],double *V[3], double a[3],double vtemp[3]) { 
        double dt,r2,r3;
        double eo =energies(&r,&V);
        r2 =*r[0]* *r[0]+*r[1]* *r[1]+*r[2]* *r[2];
        r3=r2*sqrt(r2);
        //accelaration(r,a);
        a[0] =-*r[0]/r3;
        a[1] =-*r[1]/r3;
        a[2] =-*r[2]/r3; 
              
        vtemp[0] = *V[0]+a[0]*dt*0.5;
        vtemp[1] = *V[1]+a[1]*dt*0.5;
        vtemp[2] = *V[2]+a[2]*dt*0.5;

        *r[0] += vtemp[0]*dt;
        *r[1] += vtemp[1]*dt;
        *r[2] += vtemp[2]*dt;

        //accelaration(r,a);
        a[0] =-*r[0]/r3;
        a[1] =-*r[1]/r3;
        a[2] =-*r[2]/r3; 
        *V[0] = vtemp[0]+a[0]*dt*0.5;
        *V[1] = vtemp[1]+a[1]*dt*0.5;
        *V[2] = vtemp[2]+a[2]*dt*0.5;
         printf(" %2.5e %2.5e %2.5e %2.5e \n ",r[0],r[1],r[2],(eo-energies(r,V))/eo);
        return;
    }   
    //4th degree Runge Kutta
    void rungekutta(double r[3],double V[3], double a[3],double vtemp[3] ){
        double dt,dt_param;
        double eo =energies(r,V);
        double r2 =r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
        double r3=r2*sqrt(r2);
        a[0] = -r[0]/r3;
        a[1] = -r[1]/r3;
        a[2] = -r[2]/r3; 

        //accelaration(r,a);

        double V2 =V[0]*V[0]+V[1]*V[1]+V[2]*V[2];
        dt=dt_param*sqrt(r2/V2);

        vtemp[0] = V[0]+a[0]*dt*0.5;
        vtemp[1] = V[1]+a[1]*dt*0.5;
        vtemp[2] = V[2]+a[2]*dt*0.5;

        r[0] += vtemp[0]*dt;
        r[1] += vtemp[1]*dt;
        r[2] += vtemp[2]*dt;

        //accelaration(r,a);
        a[0] = -r[0]/r3;
        a[1] = -r[1]/r3;
        a[2] = -r[2]/r3; 
        V[0] = vtemp[0]+a[0]*dt*0.5;
        V[1] = vtemp[1]+a[1]*dt*0.5;
        V[2] = vtemp[2]+a[2]*dt*0.5;
         printf(" %2.5e %2.5e %2.5e %2.5e \n ",r[0],r[1],r[2],(eo-energies(r,V))/eo);
        return ;
    }     
    int main(int argc, char **argv) {

        double r_denom,V2,r2,r3,ek,epot,cons,eo,mass,dt,tend;
        unsigned int choice;   
    double k3=0.01720209895,k2,x,y,z,i,j,k,rc,x2,y2;
    double dot1,dot2,norm_r,norm_v,norm_h,norm_N,cosE,sinE,tau;
    double eccentricity,eccentricity_vec[3],inclination,omega,w,E,M,n;
    double r[]= {-1.556577403599934E+00,-2.577116623051523E+00,-5.122853965213389E-01};    //  r_vec[3]={X,Y,Z}  V_vec[3]={Vx,Vy,Vz} ;
    double v[]= {8.012428308988644E-03,-3.979434589732452E-03,1.197182865768091E-03};
    double a[3]={0,0,0}; 
    double vtemp[3]={0,0,0} ; 
    double t=2448893.500000000,Tp=2449510.523804001976;  
    double EC= 1.869287516660393E-01 ;
    double A = 2.642546881745284E+00;
    double IN= 1.175782099005222E+01;
    double OM= 2.935665252820530E+02;
    double W = 9.747550201902160E+01;
    double MA= 2.184297687688620E+02;
    double N =  2.294404694160717E-01;

dt = 10; //DAYS
tend= 3652.5;
mass = 1.05982906e-10; // Msun
k2=k3*k3;

   norm_r = sqrt(r[0]* r[0] + r[1] * r[1] + r[2] * r[2]);
   dot1 = (r[0] * v[0]) + (r[1] * v[1]) + (r[2] * v[2]);
   dot2 = (v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]);

   for(int p=0;p<3;p++) { 
   eccentricity_vec[p] = r[p]*( dot2/k2 - (1.0/norm_r)) - v[p]* (dot1/k2); 
   }
   eccentricity  = sqrt(pow(eccentricity_vec[0],2) + pow(eccentricity_vec[1],2) + pow(eccentricity_vec[2],2));  //norm_
   printf(" e_z =%f \n",eccentricity_vec[2]);

   double semi_major = 1.0/((2.0/norm_r) - dot2/k2);

   i = r[1] * v[2] - r[2] * v[1] ;
   j = r[2] * v[0] - r[0] * v[2] ;
   k = r[0] * v[1] - r[1] * v[0];

   norm_h = sqrt(i*i+j*j+k*k);
   inclination = acos(k/norm_h)*180 /PI;

   norm_N = sqrt(j*j+i*i);
   omega = 180/PI*acos(-j/norm_N);
   w = 180/PI*acos((-j*x+i*y)/(eccentricity*norm_N));

   rc=norm_h*norm_h/k2; // rc = a*(1-eccentricity**2)
   x2 = (rc-norm_r)/eccentricity ;
   y2 = dot1/eccentricity * sqrt(rc/k2);
   cosE = (x2/semi_major)+eccentricity ;
   sinE = y2/(semi_major*sqrt(1-pow(eccentricity,2)));

   if (fabs(sinE) <= 0.707107)
   E = asin(fabs(sinE));
   if (fabs(cosE) <= 0.707107)
   E = acos(fabs(cosE));
   if (cosE >= 0 && sinE >=0) 
   E = acos(cosE);
   if (cosE<0 && sinE >= 0)
   E = PI-E ;
   if (cosE<0 && sinE < 0)
   E = PI+E ;
   if (cosE>=0 && sinE < 0) 
   E = 2*PI-E ;

   M = (E- eccentricity*sin(E))*180/PI ;
   n = 360/(pow(semi_major,1.5)*365.25); // deg/day unit
   tau = t - (M/n) ;

    printf(" eccentricity = %.16f  error = %.16f \n", eccentricity,fabs(eccentricity-EC)); 
    printf(" semi major axis = %.16f error = %.16f \n",semi_major,fabs(semi_major-A));
    printf(" inclination = %.16f error = %.16f \n",inclination,fabs(inclination-IN));
    printf(" Omega = %.16f error = %.16f \n",360-omega,fabs(OM-(360-omega)));
    printf(" omega = %.16f error = %.16f \n",w,fabs(W-w));
    printf(" E = %.16f \n",E*180/PI);
    printf(" n = %.16f error = %.16f \n",n,n-N);
    printf(" M = %.16f error = %.16f \n",M,fabs(MA-M));
    printf(" Tp= %.16f error = %.12f \n",tau,fabs(Tp-tau));


        eo =energies(r,v);
        FILE *kk;
        kk=fopen("energy5.txt","w");
        t=0.0 ;   
       // function_arr is an array of function pointers 
        void (*function_arr[])(double *,double *) = {euler,leapfrog,rungekutta}; 
        printf("enter choice : (0,1,2) \n"); 
        scanf("%d", &choice); 

        while ( t <= tend ){   // ( t <= tmax/dt )
        
        (*function_arr[choice])(r,v); 

        fprintf(kk,"%f %f %f %f %f \n",t,r[0],r[1],r[2],(eo-energies(r,v))/eo);
        printf(" %2.5e %2.5e %2.5e %2.5e \n ",r[0],r[1],r[2],(eo-energies(r,v))/eo);
        t +=dt;
        }
        fclose(kk);
        return 0; 
        }

