# dependencies

* numpy
* matplotlib

# how to use

Deze package is redelijk simpel om te gebruiken. Download eerst de TGMmodels.py file van de github paina en zet deze in de zelfde map als het bestand waar je de package in wil gebruiken. Er is een jupyter notebook op de git pagina te vinden waar alle functionaliteiten van de package in te vinden zijn. Je kan dus ook heel makkelijk de data in dit notebook vervangen met je eigen data om alle modellen op jouw data te runnen. De package bevat al een aantal differentiaal vergelijkingen maar je kan redelijk makkelijk je eigen vergelijkingen aan de package toevoegen.

Hier is een voorbeeld van een differentiaal vergelijking:
````
   class Gompertz_growth(TGmodels):
    def __init__(self, y,  vmax, c = 0.0):
        self.c = c
        self.vmax = vmax
        self.y = y
    
    def formule(self, y):
        if y <= 0:
            y = 1e-8
            self.y = 1e-8
        return self.c*y *math.log(self.vmax / y) 
    
    def params(self):
        return {"y": self.y, "c": self.c, "vmax": self.vmax}
    
    def setparams(self, y, c, vmax):
        self.c = c
        self.vmax = vmax
        self.y = y
````

Als je zelf een vergelijking toe zou willen voegen kan je een nieuwe class aan de TGMmodels.py toe voegen met de zelfde bouw als de class hierboven. Vervang dan alleen de formule en de parameters met die van de nieuwe vergelijking. Als alles goed gaat met het runnen van de voorbeeld jupyter notebook. Zou er voor elk model een BIC score uit moeten komen die aangeeft hoe goed de data op jouw model past. Hoe lager de BIC score hoe beter het model op jouw data past.

# achtergrond informatie
## Fitting methode

### Hooke & Jeeves / Direct Search

Het Hooke & Jeeves-algoritme is een eenvoudige methode om functies te optimaliseren zonder dat je afgeleiden nodig hebt. Het werkt door alleen de waarden van de functie te gebruiken.

#### Belangrijkste Eigenschappen

- **Exploratory Moves:**  
  Het algoritme past een parameter aan en kijkt of de functie beter wordt.  
  - Als de aanpassing een verbetering geeft, wordt deze geaccepteerd.  
  - Anders wordt de aanpassing teruggedraaid of kleiner gemaakt.

- **Pattern Moves:**  
  Als een goede richting is gevonden, maakt het algoritme grotere stappen in die richting om sneller een oplossing te vinden.

- **Geen afgeleiden nodig:**  
  Dit maakt het algoritme handig voor functies die moeilijk te analyseren zijn, bijvoorbeeld door ruis of plotselinge sprongen.

#### Toepassingen

- Het optimaliseren van functies waarbij de berekeningen ingewikkeld zijn.  
- Gebruikt in technische ontwerpen, machine learning (bijvoorbeeld hyperparameter-tuning), en andere experimenten.

---

## ODE-oplossers

### Runge-Kutta

De Runge-Kutta-methode is een populaire techniek om gewone differentiaalvergelijkingen op te lossen. Het biedt een goede balans tussen nauwkeurigheid en snelheid.

#### Hoe Werkt Het?

- De methode gebruikt meerdere tussenstappen om een gemiddelde helling te berekenen.  
- Deze gemiddelde helling wordt gebruikt om de oplossing bij te werken.

#### Voordelen

- **Hoge nauwkeurigheid:** Het is veel preciezer dan eenvoudigere methoden zoals Euler.  
- **Geen ingewikkelde berekeningen:** Je hoeft geen hogere-orde afgeleiden te berekenen.  
- **Flexibele stapgrootte:** De methode kan worden aangepast voor situaties waar kleinere stappen nodig zijn.

#### Toepassingen

- Veel gebruikt in wetenschappen en techniek, zoals:  
  - Simulatie van tumorgroei.  
  - Bewegingsmodellen (bijvoorbeeld planeten).  
  - Populatiemodellen.

---

### Euler-methode

De Euler-methode is de eenvoudigste manier om ODE's op te lossen. Het is minder nauwkeurig, maar makkelijk te begrijpen en te implementeren.

#### Voordelen

- **Eenvoudig:** Het is makkelijk te begrijpen en snel te implementeren.  
- **Weinig rekenkracht nodig:** Geschikt voor snelle simulaties.

#### Nadelen

- **Lage nauwkeurigheid:** Vooral bij grotere stappen is de fout groot.  
- **Kleine stapgroottes nodig:** Dit kan de berekeningen veel tijd kosten.

#### Toepassingen

- Simpele simulaties, zoals:  
  - Mechanische systemen.  
  - Basis populatiemodellen.  
- Vaak gebruikt om te leren hoe numerieke methoden werken.

---

## Begin waardes voor de modellen

Voor het kiezen van de begin waardes van de parameters van de modellen zijn een aantal dingen belangrijk.

- zet de parameters niet op 0, dit kan problemen opleveren.
- zet de parameters niet negatief, dit kan problemen opleveren.
- zet y0 nooit hoger dan ymax of lager dan ymin
- zet ymax dicht bij de hoogste waarde in je data
- zet ymin dicht bij de laagste waarde in je data
- hoe dichter bij de beste waardes je de begin waardes zet hoe sneller de package zal runnen. 