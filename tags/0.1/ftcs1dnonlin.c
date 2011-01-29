#ifdef HAVE_MALLOC_H
# include<malloc.h>
#endif
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

main()
{    FILE *fp1;
     int *iup,*idown,printed,outputinterval,lattice_size_x,i;
     float kappa,Sc,slope,den,*flux,duration,deltax,deltah,elapsedtime,max,timestep;
     float *topo,*topoold,maxchange; 

     fp1=fopen("nonlinbasedrop","w");
     lattice_size_x=102;
     deltax=1.0; /* m */
     kappa=1.0;  /* m^2/kyr */
     Sc=1.0;     /* m/m */
     iup=ivector(1,lattice_size_x);
     idown=ivector(1,lattice_size_x);
     topo=vector(1,lattice_size_x);
     topoold=vector(1,lattice_size_x);
     flux=vector(1,lattice_size_x);
     duration=1000.0; /* kyr */ 
     outputinterval=1000.0;
     maxchange=0.001;
     for (i=1;i<=lattice_size_x;i++)
      {topo[i]=1.0; 
       if ((i==lattice_size_x)||(i==lattice_size_x-1)||(i==lattice_size_x-2)) topo[i]=0.0;
       topoold[i]=topo[i];
       flux[i]=0;
       iup[i]=i+1;idown[i]=i-1;}
     iup[lattice_size_x]=lattice_size_x;idown[1]=1;
     timestep=0.0001*deltax*deltax;
     elapsedtime=0.0;printed=0;
     while (elapsedtime<duration)
      {printed=elapsedtime/outputinterval;
       for (i=1;i<=lattice_size_x;i++)
        {slope=(topoold[idown[i]]-topoold[i])/deltax;
         den=1-fabs(slope)*fabs(slope)/(Sc*Sc);
         if (den<0.01) den=0.01;
         flux[i]=kappa*slope/den;
         if (i==1) flux[i]=0.0;
         if (i==lattice_size_x) flux[i]=flux[idown[i]];
         if (flux[i]<flux[idown[i]]) flux[i]=flux[idown[i]];}
       max=0;
       for (i=1;i<lattice_size_x;i++) 
        {deltah=timestep/deltax*(flux[i]-flux[iup[i]]);
         topo[i]+=deltah;
         if (fabs(deltah)>max) max=fabs(deltah);}
       elapsedtime+=timestep;
       if (max>maxchange)
         {elapsedtime-=timestep;
          timestep/=2.0;
          for (i=1;i<lattice_size_x;i++)
           topo[i]=topoold[i];}
        else
         {if (max<0.1*maxchange) timestep*=1.2;}
       for (i=1;i<lattice_size_x;i++)
         topoold[i]=topo[i];
       if ((int)(elapsedtime/outputinterval)>printed)
        {printf("%f %f\n",elapsedtime,timestep);
         for (i=1;i<=lattice_size_x;i++)
          fprintf(fp1,"%d %f\n",i,topo[i]);
         printed=(int)(elapsedtime/outputinterval);}}
     fclose(fp1);
}  
