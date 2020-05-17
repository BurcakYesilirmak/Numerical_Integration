
    #include<stdio.h>
    #include<math.h>
    #define PI 3.1415926535897932 
    // Start=2018-09-05, Stop=2018-10-06, Step=1 d JUPITER AND SUN KEPLER PROBLEM - INTEGRATION SCHEMES

   void orbital_elements(double r_vec[3],double v_vec[3],double *e, double *a, double *inc, double *W, double *w, double *M ){
   double i,j,k,rc,x2,y2,eccentricity_vec[3],E;
   double dot1,dot2,norm_r,norm_v,norm_h,norm_N,cosE,sinE,tau;
   double t=2448893.500000000 ;  
   double k3=0.01720209895,k2 ;

   k2=k3*k3; // 1 alınırsa yörünge elemanlarının hatası çok büyük geliyor ?

    printf(" r_vec[0] = %e \n",r_vec[0]);
   norm_r = sqrt(r_vec[0]* r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]);
   dot1 = (r_vec[0] * v_vec[0]) + (r_vec[1] * v_vec[1]) + (r_vec[2] * v_vec[2]);
   dot2 = (v_vec[0] * v_vec[0]) + (v_vec[1] * v_vec[1]) + (v_vec[2] * v_vec[2]);

   for(int p=0;p<3;p++) { 
   eccentricity_vec[p] = r_vec[p]*( dot2/k2 - (1.0/norm_r)) - v_vec[p]* (dot1/k2); 
   }
   *e  = sqrt(pow(eccentricity_vec[0],2) + pow(eccentricity_vec[1],2) + pow(eccentricity_vec[2],2));  

   *a = 1.0/((2.0/norm_r) - dot2/k2);

   i = r_vec[1] * v_vec[2] - r_vec[2] * v_vec[1] ;
   j = r_vec[2] * v_vec[0] - r_vec[0] * v_vec[2] ;
   k = r_vec[0] * v_vec[1] - r_vec[1] * v_vec[0] ;

   norm_h = sqrt(i*i+j*j+k*k);
   *inc = acos(k/norm_h)*180 /PI;

   norm_N = sqrt(j*j+i*i);
   *W = acos(-j/norm_N)*180/PI;

   if (eccentricity_vec[2]<0)
   *w = 360-(acos((-j*eccentricity_vec[0]+i*eccentricity_vec[1])/(*e*norm_N))*180/PI);
   else 
   *w = acos((-j*eccentricity_vec[0]+i*eccentricity_vec[1])/(*e*norm_N))*180/PI;

   rc=norm_h*norm_h/k2; // rc = a*(1-eccentricity**2)
   x2 = (rc-norm_r)/ *e ;
   y2 = dot1/ *e * sqrt(rc/k2);
   cosE = (x2/ *a)+ *e ;
   sinE = y2/(*a*sqrt(1-pow(*e,2)));

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

   *M = (E- *e *sin(E))*180/PI ;
  // *n = 360/(pow(*a,1.5)*365.25); // deg/day unit
  // tau = t - (*M/ *n) ;
}
void position_velocity(double r_vec[3],double v_vec[3],double EC, double A, double IN, double OM, double W, double MA, double *position,double *velocity) {
    double  DEG2RAD,E,func,deriv,E_deriv;
    double X,Y,Vx,Vy,b,r;

   double k3=0.01720209895,k2 ;

   k2=k3*k3; // 1 alınırsa yörünge elemanlarının hatası çok büyük geliyor ?

    DEG2RAD =PI/180;
    IN *=DEG2RAD;
    OM *=DEG2RAD;
    W  *=DEG2RAD;
    MA *=DEG2RAD;
    double P[3],Q[3];

    P[0]= cos(OM)*cos(W)-sin(OM)*sin(W)*cos(IN);
    P[1]= sin(OM)*cos(W)+cos(OM)*sin(W)*cos(IN);
    P[2]= sin(W)*sin(IN);
    printf(" P[0]= %f P[1]= %f P[2]= %f  \n",P[0],P[1],P[2]);

    Q[0]= -sin(W)*cos(OM)-cos(W)*sin(OM)*cos(IN);
    Q[1]= -sin(W)*sin(OM)+cos(W)*cos(OM)*cos(IN);
    Q[2]= cos(W)*sin(IN);
    printf(" Q[0]= %f Q[1]= %f Q[2]= %f \n",Q[0],Q[1],Q[2]);
    
    E =MA;
    func= 1.0;
    while( fabs(func) >= 0.0000001){
    func = E - EC *sin(E) - MA ;
    deriv = 1.0 - EC*cos(E);
    E = E - func/deriv;
    }
    printf("E = %f rad \n",E);

    r = A * (1.0-EC*cos(E)) ;
    E_deriv = sqrt(k2/A)/r;
    b = A * sqrt(1.0-EC*EC);
    X = A * (cos(E)-EC);
    Y = b * sin(E);
    Vx = -A * E_deriv * sin(E);
    Vy = b * E_deriv * cos(E);

    printf("E' = %f \n",E_deriv);
    printf("r =%f X =%f Y =%f Vx =%f Vy =%f \n",r,X,Y,Vx,Vy);

    for (int i=0;i<3;i++) {
    position[i] = X * P[i] + Y * Q[i];
    velocity[i] = Vx * P[i] + Vy * Q[i];
    }
  }

 double norm(double r[3]) {   
        double norm_r = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
        return (norm_r);
    }
    double dot(double a[3],double b[3]) {
    double dot_2 = a[0]*b[0] + a[1]*b[1] + a[2]*b[2] ;
    return (dot_2);
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

    void accelaration(double r[3], double *a){
    double r2 =r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
    double r3=r2*sqrt(r2);
    for(int i=0;i<3;i++)
    a[i] = -r[i]/r3;
    return;
    } 

    void euler(double r[3],double V[3],double dt){
         double eo =energies(r,V);
         double a[3]={0,0,0};

        accelaration(r,a);
        r[0] += V[0]*dt;
        r[1] += V[1]*dt;
        r[2] += V[2]*dt;

        V[0] += a[0]*dt;   // eğer V ve r yi yerdeğiştirirsen korunum olan bir kodda
        V[1] += a[1]*dt;   // kapalı bir eğri elde edersin.
        V[2] += a[2]*dt;
        
      }

    void leapfrog(double r[3],double V[3],double dt){ 
        double t,r2,r3;
        double a[3]={0,0,0}; 
        double vtemp[3]={0,0,0} ; 
        double eo =energies(r,V);

       accelaration(r,a);
        vtemp[0] = V[0]+a[0]*dt*0.5;
        vtemp[1] = V[1]+a[1]*dt*0.5;
        vtemp[2] = V[2]+a[2]*dt*0.5;

        r[0] += vtemp[0]*dt;
        r[1] += vtemp[1]*dt;
        r[2] += vtemp[2]*dt;

        accelaration(r,a);
        V[0] = vtemp[0]+a[0]*dt*0.5;
        V[1] = vtemp[1]+a[1]*dt*0.5;
        V[2] = vtemp[2]+a[2]*dt*0.5; 
         }   

    //4th degree Runge Kutta
    void rungekutta(double r[3],double V[3],double dt_param){
        double dt,r2,V2,semi,e,inc,omega,w,M,n;
        double a[3]={0,0,0}; 
        double vtemp[3]={0,0,0} ; 
        double eo =energies(r,V);

        r2 =r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
        V2 =V[0]*V[0]+V[1]*V[1]+V[2]*V[2];
        dt=dt_param*sqrt(r2/V2);

        accelaration(r,a);
        vtemp[0] = V[0]+a[0]*dt*0.5;
        vtemp[1] = V[1]+a[1]*dt*0.5;
        vtemp[2] = V[2]+a[2]*dt*0.5;

        r[0] += vtemp[0]*dt;
        r[1] += vtemp[1]*dt;
        r[2] += vtemp[2]*dt;

        accelaration(r,a);
        V[0] = vtemp[0]+a[0]*dt*0.5;
        V[1] = vtemp[1]+a[1]*dt*0.5;
        V[2] = vtemp[2]+a[2]*dt*0.5;  
     }

    int main(int argc, char **argv) {

    double r_denom,V2,r2,r3,ek,epot,cons,eo,mass,t,dt,tend;
    unsigned int choice;   
    double e, semi, inc, omega, w, M, n;
    double position[3],velocity[3];

    //JUPITER DATA FROM NASA HORIZON
    double r[]= {-2.908569907082057E+00,-4.520452004618187E+00,8.380871524670937E-02};    //  r[3]={X,Y,Z}  V[3]={Vx,Vy,Vz} ;
    double v[]= {6.255888956995200E-03,-3.724449149570924E-03,-1.245210382528385E-04};
    double EC= 4.870777296052881E-02 ;
    double A = 5.202950141423149E+00;
    double IN= 1.303937380203285E+00;
    double OM= 1.005056413394623E+02;
    double W = 2.736878168207003E+02;
    double MA= 2.269965393247505E+02;
    double N = 8.308770556598639E-02 ;
   // double t=2448893.500000000,Tp=2459967.259820833337;  
    
    tend= 2*PI ;  
    dt= 0.0001*2*PI ;   //sqrt(size_scale*size_scale*size_scale/mass_scale); scales =1 

        // for t= 0.0 
        eo =energies(r,v);
        position_velocity(r, v, EC, A, IN, OM, W, MA,position,velocity);
        printf(" position[0]-r_vec[0] = %e \n",position[0]-r[0]);
        printf(" position[1]-r_vec[1] = %e \n",position[1]-r[1]);
        printf(" position[2]-r_vec[2] = %e \n",position[2]-r[2]);
        printf(" velocity[0]-v_vec[0] = %e \n",velocity[0]-v[0]);
        printf(" velocity[1]-v_vec[1] = %e \n",velocity[1]-v[1]);
        printf(" velocity[2]-v_vec[2] = %e \n",velocity[2]-v[2]);

        orbital_elements(r, v, &e, &semi, &inc, &omega, &w, &M);
        printf(" eccentricity = %.16e  error = %.16e \n", e,fabs(e-EC)); 
        printf(" semi major axis = %.16e error = %.16e \n",semi,fabs(semi-A));
        printf(" inclination = %.16e error = %.16e \n",inc,fabs(inc-IN));
        printf(" Omega = %.16e error = %.16e \n",omega,fabs(OM-omega));
        printf(" omega = %.16e error = %.16e \n",w,fabs(W-w));
        printf(" M = %.16e error = %.16e \n",M,fabs(MA-M));

  
        void (*function_arr[])(double *,double *,double dt) = {euler,leapfrog,rungekutta}; 
        printf("enter choice : (0,1,2) \n"); 
        scanf("%d", &choice); 
        FILE *kk;
        kk=fopen("numeric_integration.txt","w");

        t=0.0 ; 
        while ( t <= tend ){   // ( t <= tmax/dt )
        (*function_arr[choice])(r,v,dt); 
        fprintf(kk,"%e %e %e %e %e \n",t/tend,r[0],r[1],r[2],(eo-energies(r,v))/eo);
        t +=dt;
        }
        printf("  %2.5e %2.5e %2.5e %2.5e %2.5e %2.5e \n ",r[0],r[1],r[2],v[0],v[1],v[2]);
        // after 1 year
        position_velocity(r, v, EC, A, IN, OM, W, MA,position,velocity);
        printf(" position[0]-r_vec[0] = %e \n",position[0]-(-3.560786203739217E-01));
        printf(" position[1]-r_vec[1] = %e \n",position[1]-(-5.250303215848153E+00));
        printf(" position[2]-r_vec[2] = %e \n",position[2]-2.973949764781705E-02);
        printf(" velocity[0]-v_vec[0] = %e \n",velocity[0]-7.437742579081917E-03);
        printf(" velocity[1]-v_vec[1] = %e \n",velocity[1]-(-1.513445647003154E-04));
        printf(" velocity[2]-v_vec[2] = %e \n",velocity[2]-(-1.657564375930242E-04));

        orbital_elements(r, v, &e, &semi, &inc, &omega, &w, &M);
        printf(" eccentricity = %.16e  error = %.16e \n",e,fabs(e-4.861261333610650E-02)); 
        printf(" semi major axis = %.16e error = %.16e \n",semi,fabs(semi-5.204165131709652E+00 ));
        printf(" inclination = %.16e error = %.16e \n",inc,fabs(inc-1.303670294796073E+00));
        printf(" Omega = %.16e error = %.16e \n",omega,fabs(1.005163814379710E+02-omega));
        printf(" omega = %.16e error = %.16e \n",w,fabs(2.733561260937290E+02-w));
        printf(" M = %.16e error = %.16e \n",M,fabs(2.576438571877512E+02-M));

        fclose(kk);
        return (0);
    }
