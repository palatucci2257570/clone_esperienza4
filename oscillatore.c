#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M 1.0
#define EULERO 1
#define EULERO_CROMER 2
#define VELOCITY_VERLET 3
#define RK4 4
#define RK2 5

typedef struct {
  double x;
  double v;
  double a;
} Stato;
typedef struct {
  double omega2;
  double gamma;
} Parametri;
typedef double (*Accelerazione)(double t, Stato s, const Parametri *par);
typedef void (*Integratore)(double dt, Stato *s, double t, Accelerazione f, const Parametri *par);
typedef double(*posizione_esatta)(double t, Stato s0, const Parametri *par);
typedef double (*Ampiezza)(Stato s0, const Parametri *par);

double armonico(double t, Stato s, const Parametri *par){ return -(par->omega2)*(s.x);}
double energia(Stato s, const Parametri *par){return 0.5*M*((s.v*s.v)+(s.x*s.x*par->omega2));}
double errore(Stato s, double t, Stato s0, const Parametri *par, posizione_esatta x_esatto){ return fabs(x_esatto(t, s0, par)-s.x);}
double posizione_esatta_armonico(double t, Stato s0, const Parametri *par){double w = sqrt(par->omega2); return s0.x * cos((w) * t) + s0.v / (w) * sin((w) * t);}
double ampiezza_armonico(Stato s0, const Parametri *par){ return sqrt(s0.x*s0.x + s0.v*s0.v/(par->omega2)); }


void eulero(double dt, Stato *s, double t, Accelerazione f,const Parametri *par)
{
  Stato s_i=*s;
  (s->x)+=(s->v)*dt;
  (s->v)+=f(t, s_i, par)*dt;
}
void eulero_cromer(double dt, Stato *s, double t, Accelerazione f, const Parametri *par)
{
  (s->v)+=f(t, *s, par)*dt;
  (s->x)+=(s->v)*dt;
}
void velocity_verlet(double dt, Stato *s, double t, Accelerazione f, const Parametri *par)
{
  double vmid=(s->v)+f(t, *s, par)*dt*0.5; 
  (s->x)+=vmid*dt;
  (s->v)=vmid+0.5*dt*f(t, *s, par); /*se f dipende da v non va bene velocity_verlet*/
}
void rk4(double dt, Stato *s, double t, Accelerazione f, const Parametri *par)
{
  double k1[2]={(s->v), f(t, *s, par)};
  Stato s2={(s->x)+k1[0]*dt*0.5, (s->v)+0.5*dt*k1[1], 0};
  double k2[2]={s2.v, f(t+0.5*dt, s2, par)};
  Stato s3={(s->x)+k2[0]*0.5*dt, (s->v)+k2[1]*0.5*dt, 0};
  double k3[2]={s3.v, f(t+0.5*dt, s3, par)};
  Stato s4={(s->x)+k3[0]*dt, (s->v)+k3[1]*dt, 0};
  double k4[2]={s4.v, f(t+dt, s4, par)};
  (s->x)+=dt/6.0*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
  (s->v)+=dt/6.0*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
  (s->a)=f(t+dt, *s, par); /*in realtà non servirebbe aggiornarla a meno che tu non debba stampare l'intero vettore s */
}
void rk2(double dt, Stato *s, double t, Accelerazione f, const Parametri *par)
{
  double k1[2]={(s->v), f(t, *s, par)};
  Stato s2={(s->x)+k1[0]*dt*0.5, (s->v)+0.5*dt*k1[1], 0};
  double k2[2]={s2.v, f(t+0.5*dt, s2, par)};
  (s->x)+=k2[0]*dt;
  (s->v)+=k2[1]*dt;
  (s->a)=f(t+dt, *s, par);
}



int main(int argc, char *argv[])
{
  if(argc!=7){fprintf(stderr, "Uso: %s METODO DT N X0 V0 OMEGA2\n", argv[0]); return 2;}
  int metodo=atoi(argv[1]);
  double dt=atof(argv[2]);
  int n_passi=atoi(argv[3]);
  Stato s={atof(argv[4]), atof(argv[5]), 0.0};
  Stato s0=s;
  Parametri par={atof(argv[6]), 0};
  if (dt<=0 || n_passi<=0 || par.omega2<=0){fprintf(stderr, "DT N OMEGA2 DEVONO ESSERE POSITIVI\n"); return 2;}
  //scelta metodo:
  Integratore passo= NULL; 
  if (metodo == EULERO) {passo=eulero;}
  else if (metodo == EULERO_CROMER) {passo=eulero_cromer;}
  else if (metodo == VELOCITY_VERLET){ passo=velocity_verlet; }
  else if (metodo == RK4) { passo=rk4;}
  else if (metodo == RK2){ passo=rk2;}
  else {fprintf(stderr, "METODO DEVE ESSERE 1, 2, 3, 4, 5\n"); return 2;}
  //scelta tipo di sistema fisico per ora ne ho uno solo
  posizione_esatta x_esatto;
  Accelerazione a;
  Ampiezza A;
  x_esatto=posizione_esatta_armonico;
  a=armonico;
  A=ampiezza_armonico;
  
  s.a=a(0.0, s, &par);
  //stampa condizione iniziale e ciclo per calcolo e stampa stati successivi 
  fprintf(stdout, "%.10e %.10e %.10e %.10e %.10e %.10e\n", 0.,  s.x, s.v, energia(s, &par), errore(s, 0., s0, &par, x_esatto),  errore(s, 0., s0, &par, x_esatto)/A(s0, &par));
for(int n=0; n<n_passi; n++)
  {
    passo(dt, &s, n*dt, a, &par); /*qua a passo devi dare d*dt non (n+1)*dt perchè usa il tempo al passo n per calcolare il passo n+1 (le informazioni che do alla velocità/che uso per calcolare */                                              /* il passo successivo sono sempre al tempo n) */
    fprintf(stdout, "%.10e %.10e %.10e %.10e %.10e %.10e\n", (n+1)*dt, s.x, s.v, energia(s, &par), errore(s, (n+1)*dt, s0, &par, x_esatto), errore(s, (n+1)*dt, s0, &par, x_esatto)/A(s0, &par));
                                                                                                                                        //t x v E err 
                                                                                                                                        //e ci ho aggiungo l'errore diviso per l'ampiez pk mi serve nei grafici
  }
}
