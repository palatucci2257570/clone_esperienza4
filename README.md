# Esercitazione 4 - I metodi di Runge-Kutta: oscillatore armonico e pendolo smorzato

## Obiettivi

In questa esercitazione introdurremo il metodo di **Runge-Kutta del quarto ordine** (RK4), che costruiremo a partire da quello del secondo ordine (RK2). Ne verificheremo numericamente l'ordine di convergenza sull'oscillatore armonico e lo confronteremo con velocity Verlet. Riorganizzeremo inoltre il codice in modo da separare l'algoritmo di integrazione dal sistema fisico, utilizzando un **puntatore a funzione** per il calcolo dell'accelerazione. Applicheremo poi RK4 al pendolo smorzato, un sistema in cui la forza dipende dalla velocità.

Prima di cominciare fate login sul sito di GitHub.com con le vostre credenziali, poi create una cartella (ad esempio `LFC_XX`, dove `XX` è il nome del vostro gruppo), entrate al suo interno e fate il `git clone` del repository creato al momento dell'accettazione dell'esercitazione. Infine, eseguite il comando `github-inizio.sh` (vedi esercitazione 0) dal terminale che volete usare. Per assicurarvi che i *commit* siano a vostro nome, una volta entrati nella cartella clonata date i seguenti comandi:

```sh
git config user.name "Nome Cognome"
git config user.email "email-associata-a-GitHub"
```

## Consegna dell'esercitazione

La consegna deve contenere:

* `oscillatore.c`, completato a partire dallo scheletro fornito, con velocity Verlet e RK4;
* `pendolo_smorzato.c`;
* `analisi.md`, con le risposte alle domande di entrambe le parti;
* i grafici a supporto dell'analisi;
* i file `.sh` con i comandi usati per generare i dati (ad esempio `esegui_dt.sh` ed `esegui_gamma.sh`, vedi sotto);
* i file con $\Delta t$ ed $\epsilon_G$ per ciascun metodo;
* gli eventuali altri script usati per generare i dati e i grafici.

La consegna richiede un tag (vedi esercitazione 0): dopo il push che corrisponde al commit che volete consegnare, eseguite

```sh
submission_tag="submit/$(date -u +%Y%m%dT%H%M%SZ)"
git tag "$submission_tag"
git push origin "$submission_tag"
```

Il tag include solo modifiche già registrate in un commit. Per una correzione, fate un nuovo commit, push e un nuovo tag.

## 1. Runge-Kutta applicato all'oscillatore armonico

### Il sistema

Consideriamo ancora l'oscillatore armonico delle esercitazioni precedenti, con $m=1$, e la sua soluzione esatta

$$
x_{\rm esatto}(t) = x_0\cos(\omega_0 t) + \frac{v_0}{\omega_0}\sin(\omega_0 t).
$$

Come nell'esercitazione 3, misuriamo l'accuratezza di un integratore con l'errore globale sulla posizione al tempo finale $t_f$,

$$
\epsilon_G(\Delta t) = \frac{|x(t_f)-x_{\mathrm{esatto}}(t_f)|}{A},
$$

dove $A=\sqrt{x_0^2+(v_0/\omega_0)^2}$ è l'ampiezza dell'oscillazione.

### Il metodo

I metodi di Runge-Kutta si applicano a un qualunque sistema di equazioni del primo ordine nella forma

$$
\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y}),
$$

dove $\mathbf{y}$ è lo stato del sistema. Per l'oscillatore armonico $\mathbf{y}=(x, v)$ e $\mathbf{f}(t, \mathbf{y}) = (v, -\omega_0^2 x)$. A differenza di velocity Verlet, quindi, gli algoritmi non fanno distinzione tra posizioni e velocità: entrambe vengono aggiornate allo stesso modo, a partire dalle derivate calcolate da $\mathbf{f}$.

**RK2** (metodo del punto medio) stima la derivata a metà dell'intervallo e la usa per fare un passo intero. Vi servirà come primo passo verso RK4:

1. $\mathbf{k}_1 = \mathbf{f}(t_n, \mathbf{y}_n)$;
2. $\mathbf{k}_2 = \mathbf{f}\left(t_n + \frac{\Delta t}{2}, \mathbf{y}_n + \frac{\Delta t}{2}\mathbf{k}_1\right)$;
3. $\mathbf{y}_{n+1} = \mathbf{y}_n + \mathbf{k}_2\,\Delta t$.

**RK4** combina quattro stime della derivata, una all'inizio, due a metà e una alla fine dell'intervallo:

1. $\mathbf{k}_1 = \mathbf{f}(t_n, \mathbf{y}_n)$;
2. $\mathbf{k}_2 = \mathbf{f}\left(t_n + \frac{\Delta t}{2}, \mathbf{y}_n + \frac{\Delta t}{2}\mathbf{k}_1\right)$;
3. $\mathbf{k}_3 = \mathbf{f}\left(t_n + \frac{\Delta t}{2}, \mathbf{y}_n + \frac{\Delta t}{2}\mathbf{k}_2\right)$;
4. $\mathbf{k}_4 = \mathbf{f}\left(t_n + \Delta t, \mathbf{y}_n + \Delta t\,\mathbf{k}_3\right)$;
5. $\mathbf{y}_{n+1} = \mathbf{y}_n + \frac{\Delta t}{6}\left(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4\right)$.

Ogni $\mathbf{k}_i$ ha due componenti, una per $x$ e una per $v$: la prima è una velocità, la seconda un'accelerazione. I primi due coefficienti sono gli stessi in RK2 e RK4: RK4 si ottiene da RK2 aggiungendo $\mathbf{k}_3$ e $\mathbf{k}_4$ e cambiando l'ultimo passo.

### Analisi e consegne

Rispondete alle seguenti domande nella prima parte del file `analisi.md`:

1.1 Quale dipendenza dell'errore $\epsilon_G$ da $\Delta t$ osservate per velocity Verlet e RK4? Qual è la pendenza dei dati in un grafico log-log?

1.2 Quale algoritmo è più indicato per studiare l'oscillatore armonico su tempi lunghi: velocity Verlet o RK4? Integrate l'oscillatore per molti periodi con entrambi i metodi e graficate la variazione relativa dell'energia $(E(t)-E(0))/E(0)$ in funzione del tempo. Come si comporta l'energia nei due casi? A parità di $\Delta t$ l'errore $\epsilon_G$ di RK4 è molto più piccolo di quello di velocity Verlet: esiste un tempo oltre il quale RK4 conserva l'energia peggio di velocity Verlet?

Supportate le vostre considerazioni con il grafico log-log di $\epsilon_G$ e con il grafico dell'energia in funzione del tempo, da allegare alla consegna.

### Estensioni opzionali

Si consiglia di passare al punto 2 dopo non più di un'ora e mezzo dall'inizio dell'esercitazione, anche se i punti opzionali non sono stati completati.

1.3 Per passi di integrazione molto piccoli, l'errore continua indefinitamente a seguire la stessa legge? Se no, come interpretate quello che osservate?

1.4 Un passo di RK4 richiede quattro valutazioni dell'accelerazione, uno di velocity Verlet una sola. Rifate il grafico di $\epsilon_G$ in funzione del numero totale di valutazioni della forza invece che di $\Delta t$: a parità di costo, quale metodo è più accurato?

### Implementazione del codice

Nel repository trovate lo scheletro `oscillatore.c`, da completare seguendo i commenti `TASK` e `TODO`. **Interfaccia e output del programma restano quelli delle esercitazioni 2 e 3**:

```text
./oscillatore METODO DT N X0 V0 OMEGA2
```

con `OMEGA2` $=\omega_0^2$ e, su stdout, $N+1$ righe (la prima per la condizione iniziale) con le cinque colonne

```text
t x v E errore
```

dove `errore` è $|x - x_{\rm esatto}(t)|$, come nell'esercitazione 3. I valori di `METODO` restano quelli dell'esercitazione 3 per velocity Verlet (3), a cui si aggiunge RK4 (4). Nello scheletro è già presente Eulero (1), come esempio. Come nelle esercitazioni precedenti, il programma deve terminare con un messaggio di errore e un codice diverso da 0 se il numero di argomenti è sbagliato, se `METODO` non corrisponde a un metodo implementato, o se `DT`, `N` o `OMEGA2` non sono positivi. Per compilare usate `make`.

L'elemento nuovo è la separazione tra integratore e sistema fisico. Gli integratori non devono contenere nulla di specifico dell'oscillatore: tutta la fisica va in una funzione che, dati il tempo e lo stato del sistema, restituisce l'accelerazione. L'integratore riceve questa funzione tramite un **puntatore a funzione**. Lo scheletro è organizzato così:

```c
typedef struct { /* ... */ } Stato;      /* posizione x, velocità v, accelerazione a */
typedef struct { /* ... */ } Parametri;  /* parametri del sistema fisico */

/* puntatore a una funzione che, dati t e lo stato s, restituisce l'accelerazione */
typedef double (*Accelerazione)(double t, Stato s, const Parametri *par);

/* puntatore a un integratore: avanza *s di un passo dt a partire dal tempo t */
typedef void (*Integratore)(double t, double dt, Stato *s, const Parametri *par, Accelerazione f);

void eulero(double t, double dt, Stato *s, const Parametri *par, Accelerazione f);  /* già scritto */
void velocity_verlet(double t, double dt, Stato *s, const Parametri *par, Accelerazione f);
void rk4(double t, double dt, Stato *s, const Parametri *par, Accelerazione f);
```

I passi da seguire sono:

1. **TASK 1.** Completate le strutture `Stato`, con i campi `x`, `v` e `a` (l'accelerazione nello stato corrente), e `Parametri`, con i parametri del sistema (per l'oscillatore $\omega_0^2$; per la seconda parte servirà anche $\gamma$). Scrivete poi la funzione `armonico`, che ha la forma richiesta dal tipo `Accelerazione` e restituisce $a = -\omega_0^2 x$, e l'energia.
2. **TASK 2.** Il metodo di Eulero è già scritto nella nuova forma, come esempio: l'integratore non conosce il sistema e chiede l'accelerazione a `f`. Riscrivete allo stesso modo velocity Verlet a partire dal vostro codice dell'esercitazione 3. Il grande vantaggio di velocity Verlet è che richiede **una sola valutazione dell'accelerazione per passo**: l'accelerazione nella vecchia posizione è quella calcolata al passo precedente, che trovate in `s->a`. Chiamate `f` una sola volta, nella nuova posizione, e salvate il risultato in `s->a` per il passo successivo. Per questo, nel `main`, l'accelerazione iniziale va calcolata prima del ciclo.
3. **TASK 3.** Scrivete `rk4` in due fasi. Scrivete prima RK2 ($\mathbf{k}_1$, $\mathbf{k}_2$ e il passo con $\mathbf{k}_2$) e controllate che funzioni: ad esempio, che l'errore $\epsilon_G$ scali come $\Delta t^2$. Poi aggiungete $\mathbf{k}_3$ e $\mathbf{k}_4$ e cambiate l'ultimo passo, trasformandolo in RK4. Per ciascun $\mathbf{k}_i$ la componente per $x$ è una velocità e quella per $v$ un'accelerazione, calcolata con `f`.
4. **TASK 4.** Nel `main`, scegliete l'integratore **una sola volta**, prima del ciclo, assegnando al puntatore `passo` (di tipo `Integratore`) la funzione corrispondente a `METODO`. Nel ciclo chiamate poi sempre `passo(...)`, passando `armonico` come ultimo argomento: il ciclo non deve sapere quale metodo state usando.

Con questa struttura, per integrare un sistema diverso basta scrivere una nuova funzione per l'accelerazione, senza toccare gli integratori. È quello che farete nella seconda parte.

**Consiglio:** prima di studiare l'errore, controllate che velocity Verlet dia esattamente lo stesso output dell'esercitazione 3, e che RK4 funzioni: ad esempio, guardate la traiettoria nello spazio delle fasi e verificate che l'energia sia (quasi) costante.

### Generazione dell'output

Come nell'esercitazione 3, per studiare l'errore in funzione di $\Delta t$ dovete eseguire il programma più volte, una per ogni metodo e per ogni $\Delta t$, salvando ogni volta l'output in un file diverso. Il programma esegue una sola simulazione: non aggiungete cicli su $\Delta t$ o sui metodi dentro al codice C. Riportate invece i comandi in un file di testo con estensione `.sh`, ad esempio `esegui_dt.sh`, con una riga per ogni esecuzione:

```sh
./oscillatore 3 0.001 10240 1 0.5 4 > vv_dt0.001.dat
./oscillatore 3 0.002 5120 1 0.5 4 > vv_dt0.002.dat
...
./oscillatore 4 0.001 10240 1 0.5 4 > rk4_dt0.001.dat
./oscillatore 4 0.002 5120 1 0.5 4 > rk4_dt0.002.dat
...
```

e lanciatelo con

```sh
bash esegui_dt.sh
```

**Consigli:**

1. Come nell'esercitazione 3, potete usare $\Delta t$ = 0.001, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, con $N = 10240$ per $\Delta t = 0.001$ e dimezzando $N$ ogni volta che raddoppiate $\Delta t$, in modo che il tempo finale $t_f = N\,\Delta t = 10.24$ sia lo stesso per tutte le esecuzioni.
2. Per la domanda 1.3 aggiungete valori di $\Delta t$ più piccoli (e $N$ più grandi), ad esempio fino a $\Delta t = 10^{-4}$: l'errore diventa molto piccolo, soprattutto con RK4, quindi stampate abbastanza cifre (ad esempio con `%.10e`).

### Come rispondere alle domande

Valgono gli stessi suggerimenti dell'esercitazione 3: da ciascun file prendete la colonna `errore` dell'**ultima riga**, cioè l'errore al tempo finale $t_f$, e dividetela per l'ampiezza $A$: ottenete $\epsilon_G(\Delta t)$. Raccogliete le coppie $\Delta t$, $\epsilon_G$ in un file per ciascun metodo (ad esempio `errore_rk4.dat`) e graficatele in scala log-log. Stimate la pendenza con un fit lineare dei logaritmi, oppure confrontando i dati con rette di pendenza nota. Rappresentate tutti i metodi sullo stesso grafico.

Per la domanda 1.4, il numero totale di valutazioni dell'accelerazione è $N$ per velocity Verlet (più quella iniziale) e $4N$ per RK4.

**Energia su tempi lunghi (domanda 1.2).** Usate un $\Delta t$ non troppo piccolo rispetto al periodo, altrimenti per vedere le differenze servono simulazioni lunghissime. Ad esempio, con $\omega_0 = 1$ (periodo $T_0 = 2\pi$), $\Delta t = 0.2$ e $N = 20000$ si integrano circa 640 periodi:

```sh
./oscillatore 3 0.2 20000 1 0 1 > vv_lungo.dat
./oscillatore 4 0.2 20000 1 0 1 > rk4_lungo.dat
```

(anche questi comandi vanno in un file `.sh`). Distinguete tra un errore sull'energia che oscilla restando limitato e uno che cresce, o diminuisce, sistematicamente nel tempo. Il comportamento di velocity Verlet è legato al fatto che è un integratore **simplettico**: cercate sul libro o sulle note del corso che cosa significa. Attenzione: la conservazione dell'energia non dice tutto. Confrontate anche, ad esempio, la colonna `errore` alla fine della simulazione.

## 2. Il pendolo smorzato

### Il sistema

Consideriamo un pendolo di lunghezza $L$ e massa $m$ soggetto, oltre alla gravità, a una forza dissipativa proporzionale alla velocità, come succede per un oggetto che si muove lentamente in un fluido viscoso. In forma di sistema del primo ordine, il moto è descritto da

$$
\begin{cases}
\frac{d\theta}{dt}=\omega\\
\frac{d\omega}{dt}=-\frac{g}{L}\sin\theta - \gamma\,\omega,
\end{cases}
$$

dove $\gamma \geq 0$ è il coefficiente di smorzamento, che ha le dimensioni dell'inverso di un tempo. L'energia meccanica è, come in precedenza,

$$
E= \frac{1}{2} m L ^2\omega^2 + mg L (1-\cos\theta),
$$

ma per $\gamma > 0$ non è più una quantità conservata.

Nel limite delle piccole oscillazioni, $\sin\theta \simeq \theta$, il sistema diventa un oscillatore armonico smorzato con pulsazione propria $\omega_0=\sqrt{g/L}$. In questo limite si distinguono tre regimi, a seconda del valore di $\gamma$ rispetto a $2\omega_0$: smorzamento sottocritico ($\gamma < 2\omega_0$), critico ($\gamma = 2\omega_0$) e sovracritico ($\gamma > 2\omega_0$). Trovate la soluzione esatta nei tre casi nel libro o sulle note del corso.

### Il metodo

Useremo RK4, con $\mathbf{y} = (\theta, \omega)$ e

$$
\mathbf{f}(t, \mathbf{y}) = \left(\omega, -\frac{g}{L}\sin\theta - \gamma\,\omega\right).
$$

Grazie alla separazione tra integratore e sistema fisico introdotta nella prima parte, per passare dall'oscillatore al pendolo smorzato basta cambiare la funzione che calcola l'accelerazione.

### Analisi e consegne

Rispondete alle seguenti domande nella seconda parte del file `analisi.md`:

2.1 Come cambiano il moto $\theta(t)$ e la traiettoria nello spazio delle fasi $(\theta, \omega)$ all'aumentare di $\gamma$? Confrontate i risultati con il caso conservativo $\gamma = 0$.

2.2 Come evolve l'energia meccanica nel tempo? Come dipende il suo andamento da $\gamma$?

2.3 Perché in questa parte abbiamo usato Runge-Kutta e non velocity Verlet? Per rispondere, nel vostro `pendolo_smorzato.c` sostituite la chiamata a `rk4` con una chiamata a `velocity_verlet`, copiato senza modifiche dalla prima parte (ricordatevi di calcolare `s.a` prima del ciclo): grazie al puntatore a funzione basta cambiare il nome dell'integratore. Confrontate i risultati con quelli di Runge-Kutta per $\gamma = 0$ e per $\gamma > 0$. C'è un motivo per preferire velocity Verlet nel caso dissipativo? Ricordatevi di ripristinare RK4 prima della consegna.

Supportate le vostre considerazioni con dei grafici (almeno $\theta(t)$, lo spazio delle fasi ed $E(t)$ per diversi valori di $\gamma$, compreso $\gamma = 0$).

### Estensioni opzionali

2.4 **Smorzamento critico.** Partendo da un angolo iniziale piccolo e da fermo, cercate numericamente il valore di $\gamma$ per cui il pendolo torna più rapidamente nella posizione di equilibrio senza oltrepassarla, e confrontatelo con $2\omega_0$. Cosa cambia se l'angolo iniziale è grande?

2.5 **Bilancio energetico.** Dimostrate a partire dalle equazioni del moto che $dE/dt = -\gamma m L^2\omega^2$. Verificate numericamente questa relazione, confrontando $E(t)$ con $E(0) - \int_0^t \gamma m L^2\omega^2(t')\,dt'$, dove l'integrale può essere calcolato, ad esempio, con la regola dei trapezi.

### Implementazione del codice

Scrivete un programma C, il cui file sorgente **deve chiamarsi** `pendolo_smorzato.c`, che integri la dinamica mediante RK4. Potete partire dal vostro `pendolo.c` dell'esercitazione 3. I tipi `Stato`, `Parametri` e `Accelerazione` e la funzione `rk4` devono essere quelli della prima parte, copiati senza modifiche: se avete separato correttamente integratore e sistema fisico, dovrete scrivere solo la funzione da integrare,

```c
double pendolo(double t, Stato s, const Parametri *par);  /* -(g/L) sin(theta) - gamma omega */
```

che include lo smorzamento: con `gamma` = 0 si ritrova il pendolo non smorzato dell'esercitazione 3.

In questo caso il campo `x` di `Stato` contiene $\theta$, il campo `v` contiene $\omega$ e il campo `a` l'accelerazione angolare, mentre `Parametri` contiene $g/L$, al posto di $\omega_0^2$, e il coefficiente di smorzamento $\gamma$.

Il programma deve essere eseguito con la sintassi

```text
./pendolo_smorzato dt nsteps stride theta0 omega0 length gamma
```

dove gli argomenti sono gli stessi di `pendolo.c` dell'esercitazione 3, con in aggiunta `gamma`, che è il coefficiente di smorzamento $\gamma$, in $\mathrm{s^{-1}}$. Il programma deve verificare che sia stato fornito il numero corretto di argomenti e, in caso contrario, terminare con un messaggio di errore e un codice diverso da 0.

Come nell'esercitazione 3, fissate $m=1$ e definite $g = 9.81\ \mathrm{m/s^2}$ come costante nel codice. Per compilare usate `make` (il `Makefile` compila sia `oscillatore` sia `pendolo_smorzato`).

### Generazione dell'output

Valgono le stesse convenzioni del pendolo dell'esercitazione 3: va stampata sempre la condizione iniziale ($t=0$) e poi una riga ogni `stride` passi. Ogni riga deve contenere esattamente quattro colonne separate da spazi, nell'ordine

```text
t theta omega E
```

senza intestazioni o altre informazioni.

Anche qui ogni esecuzione corrisponde a un solo valore di $\gamma$ (e di $\theta_0$, $\Delta t$, ...). Per esplorare più valori raccogliete i comandi in un file `.sh`, ad esempio `esegui_gamma.sh`:

```sh
./pendolo_smorzato 0.001 20000 10 1 0 1 0    > gamma0.dat
./pendolo_smorzato 0.001 20000 10 1 0 1 0.16 > gamma0.16.dat
./pendolo_smorzato 0.001 20000 10 1 0 1 0.63 > gamma0.63.dat
...
```

e lanciatelo con `bash esegui_gamma.sh`.

### Come rispondere alle domande

**Scelta dei parametri.** Conviene esprimere $\gamma$ in unità della pulsazione propria $\omega_0 = \sqrt{g/L}$ (per $L = 1$ m, $\omega_0 \simeq 3.13\ \mathrm{s^{-1}}$), scegliendo valori sia molto più piccoli sia più grandi di $2\omega_0$: ad esempio $\gamma/\omega_0 = 0, 0.05, 0.2, 1, 2, 5$. Come nell'esercitazione 3, il passo di integrazione va scelto piccolo rispetto alle scale di tempo del sistema, cioè $T_0 = 2\pi/\omega_0$ e $1/\gamma$. Integrate per un tempo sufficiente a vedere il moto esaurirsi, almeno per i valori di $\gamma$ più grandi.

**Moto e spazio delle fasi.** Per confrontare i diversi regimi graficate $\theta(t)$ per diversi valori di $\gamma$ sullo stesso grafico, e fate lo stesso per le traiettorie nello spazio delle fasi. Chiedetevi verso quale punto si dirigono le traiettorie e se il pendolo oscilla attorno alla posizione di equilibrio prima di raggiungerla.

**Energia.** Graficate $E(t)$ per diversi valori di $\gamma$. Per capire con quale legge decresce l'energia può essere utile una scala semilogaritmica (logaritmica solo sull'asse $y$), in cui un decadimento esponenziale appare come una retta. Osservate con attenzione anche l'andamento su tempi più brevi di un periodo: l'energia diminuisce con la stessa rapidità in tutti gli istanti?

**Velocity Verlet con lo smorzamento (domanda 2.3).** Velocity Verlet calcola l'accelerazione nella nuova posizione, $a_{n+1}$, prima di conoscere la nuova velocità $v_{n+1}$: chiedetevi con quale velocità viene calcolata, nel vostro codice, quando l'accelerazione dipende da $\omega$. Per confrontare i due integratori, ripetete la stessa simulazione con due valori di $\Delta t$ (ad esempio $10^{-3}$ e $2\cdot 10^{-3}$) e guardate di quanto cambia $\theta$ a un tempo fissato: per un metodo di ordine $p$ la differenza scala come $\Delta t^p$.

## Per chi ha finito: uno script bash con un ciclo

Se avete completato tutta l'esercitazione, comprese le estensioni opzionali, e non l'avete già fatto nell'esercitazione 3, cercate su internet che cos'è uno *script bash* e come si scrive un ciclo `for` in bash. Provate poi a riscrivere `esegui_dt.sh` con due cicli annidati, uno sui metodi e uno sui valori di $\Delta t$, che a ogni iterazione calcolano $N$, lanciano `./oscillatore` e salvano l'output in un file con un nome diverso. Controllate che i file prodotti siano identici a quelli ottenuti con la versione precedente.
