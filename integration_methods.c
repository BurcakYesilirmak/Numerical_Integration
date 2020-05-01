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
        double norm_r,r_denom,V2,rx,ry,rz,Vx,Vy,Vz,r2,r3;
        double t,dt=0.01,ek,epot,cons,eo;
        int tend=100;
        unsigned int choice;

        printf("( Enter components of (rx,ry,rz) : \n");
        scanf("%lf %lf %lf",&rx,&ry,&rz);
        printf("( Enter components of (Vx,Vy,Vz) : \n");
        scanf("%lf %lf %lf",&Vx,&Vy,&Vz);
        printf(" Enter time interval :"); 
        scanf("%lf",&dt); //0.01 ile başla
        double r[]= {rx,ry,rx};  
        double V[]= {Vx,Vy,Vz};  
        double a[]={0,0,0}; 
        double vtemp[]={0,0,0} ; 

        eo =energies(r,V);


        FILE *kk;
        kk=fopen("energy5.txt","w");
        t=0.0 ;   
       // function_arr is an array of function pointers 
        void (*function_arr[])(double *,double *,double *,double *) = {euler,leapfrog,rungekutta}; 
        printf("enter choice : (0,1,2) \n"); 
        scanf("%d", &choice); 

        while ( t <= tend ){   // ( t <= tmax/dt )
        
        (*function_arr[choice])(&r,&V,a,vtemp); 

        fprintf(kk,"%f %f %f %f %f \n",t,r[0],r[1],r[2],(eo-energies(r,V))/eo);
        t +=dt;
        }
        fclose(kk);
        return 0; 
        }
