#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M 1.0
#define G 9.805
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
  double omega2; // stiamo considerando omega2 come omega0 (non al quadrato) del mio oscillatore
  double gamma;
} Parametri;
typedef double (*Accelerazione)(double t, Stato s, const Parametri *par);
typedef void (*Integratore)(double dt, Stato *s, double t, Accelerazione f, const Parametri *par);
typedef double(*posizione_esatta)(double t, Stato s0, const Parametri *par);
typedef double (*Ampiezza)(Stato s0, const Parametri *par);

double smorzato(double t, Stato s, const Parametri *par){ return (-(par -> omega2)*(par -> omega2)*sin(s.x)-(par -> gamma)*(s.v)) ;} // abbiamo inserito la funzione di tipo accelerazione "smorzato"
double energia(Stato s, const Parametri *par){return M*(0.5*(G*G*(s.v*s.v)/((par->omega2)*(par->omega2)*(par->omega2)*(par->omega2))+(G*G/((par->omega2)*(par->omega2)*(1-cos(s.x))))));}
double errore(Stato s, double t, Stato s0, const Parametri *par, posizione_esatta x_esatto){ return fabs(x_esatto(t, s0, par)-s.x);}
double posizione_esatta_armonico(double t, Stato s0, const Parametri *par){double w = sqrt(par->omega2); return s0.x * cos((w) * t) + s0.v / (w) * sin((w) * t);}
double ampiezza_(Stato s0, const Parametri *par){ return sqrt(s0.x*s0.x + s0.v*s0.v/(par->omega2)); }



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
  if(argc!=8){fprintf(stderr, "Uso: %s DT N STRIDE THETA0 OMEGA0 LENGHT GAMMA\n", argv[0]); return 2;}
  double dt=atof(argv[1]);
  int n_passi=atoi(argv[2]);
  int stride=atoi(argv[3]);
  Stato s={atof(argv[4]), atof(argv[5]), 0.0};
  if (dt<=0 || n_passi<=0 || atof(argv[6]) <=0 || atoi(argv[3])<=0){fprintf(stderr, "DT N STRIDE LENGHT DEVONO ESSERE POSITIVI\n"); return 2;}
  Parametri par={sqrt(G/atof(argv[6])), atof(argv[7])}; // abbiamo messo i parametri richiesti
  //scelta metodo:
    Integratore passo= rk4; // abbiamo messo un solo integratore

  //scelta tipo di sistema fisico per ora ne ho uno solo
 
  s.a=smorzato(0.0, s, &par);
  //stampa condizione iniziale e ciclo per calcolo e stampa stati successivi 
  fprintf(stdout, "%.10e %.10e %.10e %.10e\n", 0.,  s.x, s.v, energia(s, &par));
for(int n=0; n<n_passi; n++)
  {
    passo(dt, &s, n*dt, smorzato, &par); /*qua a passo devi dare d*dt non (n+1)*dt perchè usa il tempo al passo n per calcolare il passo n+1 (le informazioni che do alla velocità/che uso per calcolare */                                              /* il passo successivo sono sempre al tempo n) */
    if (n % stride == 0) {
      fprintf(stdout, "%.10e %.10e %.10e %.10e\n", (n+1)*dt, s.x, s.v, energia(s, &par));
    }
                                                                                                                                        //t x v E err 
                                                                                                                                        //e ci ho aggiungo l'errore diviso per l'ampiez pk mi serve nei grafici
  }
}

