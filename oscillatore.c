#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define EULERO 1
#define VELOCITY_VERLET 3
#define RK4 4
#define M 1.0

/* TASK 1: lo stato del sistema, con i campi x (posizione), v (velocita') e
 * a (accelerazione nello stato corrente). Il campo a serve a velocity Verlet
 * per riusare al passo successivo l'accelerazione gia' calcolata. */
typedef struct {
  double x;
  double v;
  double a;
} Stato;

/* TASK 1: i parametri del sistema fisico, con i campi omega2 e gamma. Per
 * l'oscillatore serve solo omega2; gamma servira' per il pendolo smorzato
 * della seconda parte (qui vale 0). */
typedef struct {
  double omega2;
  double gamma;
} Parametri;

/* Accelerazione e' il nome di un tipo: "puntatore a una funzione che riceve
 * (double t, Stato s, const Parametri *par) e restituisce un double",
 * l'accelerazione dv/dt nello stato s. Una variabile di questo tipo puo'
 * contenere l'indirizzo di funzioni come armonico (vedi sotto): gli
 * integratori la ricevono come argomento f e la chiamano con f(t, s, par)
 * senza sapere quale sistema stanno integrando. */
typedef double (*Accelerazione)(double t, Stato s, const Parametri *par); 

/* Puntatore a un integratore: una funzione che avanza lo stato *s di un passo
 * dt, a partire dal tempo t, usando solo la funzione f per conoscere
 * l'accelerazione. Tutti gli integratori qui sotto hanno questa forma. */
typedef void (*Integratore)(double t, double dt, Stato *s, const Parametri *par, Accelerazione f); //////////////////////////////////////////////////////////////////////////////////////

/* TASK 1: l'accelerazione dell'oscillatore armonico, a = -omega^2 x. Ha la
 * forma richiesta dal tipo Accelerazione. Sostituite il valore provvisorio. */
double armonico(double t, Stato s, const Parametri *par)
{
  return -(par->omega2)*(s.x); // (*par).omega2
}

/* TASK 1: energia meccanica, E = m v^2 / 2 + m omega^2 x^2 / 2.
 * Potete riusare la funzione dell'esercitazione 3. */
double energia(Stato s, const Parametri *par)
{
  return 0.5*(M*(s.v)*(s.v)+M*(par->omega2)*(s.x)*(s.x)); // TODO
}

/* Esempio: un passo del metodo di Eulero scritto nella nuova forma.
 * L'integratore non sa quale sistema sta integrando: chiede
 * l'accelerazione a f. */
void eulero(double t, double dt, Stato *s, const Parametri *par, Accelerazione f)
{
    double a = f(t, *s, par);
    s->x = s->x + s->v * dt;
    s->v = s->v + a * dt;
}

/* TASK 2: un passo di velocity Verlet nella nuova forma, partendo dal vostro
 * codice dell'esercitazione 3. L'accelerazione nella vecchia posizione e' gia'
 * in s->a: chiamate f una sola volta, nella nuova posizione, e salvate il
 * risultato in s->a per il passo successivo. */
void velocity_verlet(double t, double dt, Stato *s, const Parametri *par, Accelerazione f)
{
  (s->x)+=(s->v)* dt;
  double vmetà=(s->v)+0.5*(s->a)*dt;
  (s->v)=vmetà+0.5*dt*f(t, *s, par);
}

/* TASK 3: un passo di Runge-Kutta. Non deve contenere nulla di specifico
 * dell'oscillatore: solo chiamate a f. Procedete in due fasi:
 *  a) scrivete prima RK2 (k1 e k2, poi il passo con k2) e controllate che
 *     funzioni, ad esempio che l'errore scali come dt^2;
 *  b) aggiungete k3 e k4 e cambiate l'ultimo passo, trasformandolo in RK4.
 * Ogni k ha due componenti: quella per x (una velocita') e quella per v
 * (un'accelerazione, calcolata con f). */
void rk4(double t, double dt, Stato *s, const Parametri *par, Accelerazione f)
{
  double k1[2]={s->v, f(t, *s, par)};
  k2_x=k1[0]+k1[1]*dt*0.5;
  double xmid=(s->x)+0.5*k2_x*dt;
  k2_v=f(t+0.5*dt, 
  double k2[2]={k1[0]+k1[1]*dt*0.5, f(t+0.5*dt, *s, *par)};
  (s->x)+=k2[0]*dt;
  (s->v)+=k2[1]*dt;
}

int main(int argc, char *argv[])
{
    /* Il programma si usa come nelle esercitazioni 2 e 3:
     *     ./oscillatore METODO DT N X0 V0 OMEGA2
     * Potete sostituire questa parte con il vostro codice dell'esercitazione 3. */
    if (argc != 7) {
        fprintf(stderr, "Uso: %s METODO DT N X0 V0 OMEGA2\n", argv[0]);
        return 2;
    }

    int metodo = atoi(argv[1]);
    double dt = atof(argv[2]);
    int n_passi = atoi(argv[3]);
    double x0 = atof(argv[4]);
    double v0 = atof(argv[5]);
    Parametri par = {atof(argv[6]), 0.0};
    if (dt <= 0.0 || n_passi <= 0 || par.omega2 <= 0.0) {
        fprintf(stderr, "Errore: DT, N e OMEGA2 devono essere positivi\n");
        return 2;
    }

    /* TASK 4: scegliete l'integratore una volta sola, prima del ciclo,
     * assegnando al puntatore a funzione `passo` la funzione giusta. */
    Integratore passo = NULL;
    if (metodo == EULERO) {
        passo = eulero;
    }
    // TODO: aggiungete un else if per VELOCITY_VERLET e uno per RK4
    else {
        fprintf(stderr, "Errore: METODO deve essere 1, 3 o 4\n");
        return 2;
    }

    /* Lo stato iniziale. L'accelerazione iniziale va calcolata una volta,
     * prima del ciclo: velocity Verlet la usa al primo passo. */
    Stato s = {x0, v0, 0.0};
    s.a = armonico(0.0, s, &par);

    /* TASK 4: l'output resta quello delle esercitazioni 2 e 3: N + 1 righe
     * con le cinque colonne
     *     t x v E errore
     * per esempio con
     *     printf("%.10e %.10e %.10e %.10e %.10e\n", t, s.x, s.v, E, errore);
     * Stampate prima la condizione iniziale (t = 0), poi nel ciclo fate N
     * passi chiamando l'integratore tramite il puntatore,
     *     passo(n * dt, dt, &s, &par, armonico);
     * e stampate lo stato al nuovo tempo t = (n + 1) * dt.
     * Il ciclo non deve sapere quale metodo e' stato scelto. */
    // TODO: stampa della condizione iniziale
    for (int n = 0; n < n_passi; ++n) {
        // TODO
    }

    return 0;
}
