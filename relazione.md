1.1) rk4 ha pendenza 3, velocity-verlet 2. Si allega il grafico:
![dipendeza dell'errore normalizzato da dt](errore_confronto.png)
1.2) Come si può vedere dal grafico, per periodi di oscillazioni molto (ma davvero molto) lunghi (dichiarati nel bash allegato), l'energia di RK4 si conserva peggio di Velocity Verlet: ciò è dovuto al fatto che RK4 non è simplettico. 
![energia vs tempo](energia_vs_t.png) 
