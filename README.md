# Tumor growth models
## achtergrond informatie
### Fitting methode

#### Hooke & Jeeves / Direct Search

Het Hooke & Jeeves-algoritme is een eenvoudige methode om functies te optimaliseren zonder dat je afgeleiden nodig hebt. Het werkt door alleen de waarden van de functie te gebruiken.

##### Belangrijkste Eigenschappen

- **Exploratory Moves:**  
  Het algoritme past een parameter aan en kijkt of de functie beter wordt.  
  - Als de aanpassing een verbetering geeft, wordt deze geaccepteerd.  
  - Anders wordt de aanpassing teruggedraaid of kleiner gemaakt.

- **Pattern Moves:**  
  Als een goede richting is gevonden, maakt het algoritme grotere stappen in die richting om sneller een oplossing te vinden.

- **Geen afgeleiden nodig:**  
  Dit maakt het algoritme handig voor functies die moeilijk te analyseren zijn, bijvoorbeeld door ruis of plotselinge sprongen.

##### Toepassingen

- Het optimaliseren van functies waarbij de berekeningen ingewikkeld zijn.  
- Gebruikt in technische ontwerpen, machine learning (bijvoorbeeld hyperparameter-tuning), en andere experimenten.

---

### ODE-oplossers

#### Runge-Kutta

De Runge-Kutta-methode is een populaire techniek om gewone differentiaalvergelijkingen op te lossen. Het biedt een goede balans tussen nauwkeurigheid en snelheid.

##### Hoe Werkt Het?

- De methode gebruikt meerdere tussenstappen om een gemiddelde helling te berekenen.  
- Deze gemiddelde helling wordt gebruikt om de oplossing bij te werken.

##### Voordelen

- **Hoge nauwkeurigheid:** Het is veel preciezer dan eenvoudigere methoden zoals Euler.  
- **Geen ingewikkelde berekeningen:** Je hoeft geen hogere-orde afgeleiden te berekenen.  
- **Flexibele stapgrootte:** De methode kan worden aangepast voor situaties waar kleinere stappen nodig zijn.

##### Toepassingen

- Veel gebruikt in wetenschappen en techniek, zoals:  
  - Simulatie van tumorgroei.  
  - Bewegingsmodellen (bijvoorbeeld planeten).  
  - Populatiemodellen.

---

#### Euler-methode

De Euler-methode is de eenvoudigste manier om ODE's op te lossen. Het is minder nauwkeurig, maar makkelijk te begrijpen en te implementeren.

##### Voordelen

- **Eenvoudig:** Het is makkelijk te begrijpen en snel te implementeren.  
- **Weinig rekenkracht nodig:** Geschikt voor snelle simulaties.

##### Nadelen

- **Lage nauwkeurigheid:** Vooral bij grotere stappen is de fout groot.  
- **Kleine stapgroottes nodig:** Dit kan de berekeningen veel tijd kosten.

##### Toepassingen

- Simpele simulaties, zoals:  
  - Mechanische systemen.  
  - Basis populatiemodellen.  
- Vaak gebruikt om te leren hoe numerieke methoden werken.

---

### Begin waardes voor de modellen

Voor het kiezen van de begin waardes van de parameters van de modellen zijn een aantal dingen belangrijk.

- zet de parameters niet op 0, dit kan problemen opleveren.
- zet de parameters niet negatief, dit kan problemen opleveren.
- zet y0 nooit hoger dan ymax of lager dan ymin
- zet ymax dicht bij de hoogste waarde in je data
- zet ymin dicht bij de laagste waarde in je data
- hoe dichter bij de beste waardes je de begin waardes zet hoe sneller de package zal runnen. 
## User's guide   
  
### Step 1   
Om te starten moet het model worden gedefineerd, dit kan van de volgende:  
  
``
import TGModels as model  
    [  
        model.Lineair_growth,  
        model.Exponentieel_toenemende_groei,  
        model.Exponentieel_afvlakkende_groei,  
        model.Mendelsohn_growth,  
        model.Allee_effect_growth,  
        model.Lineair_gelimiteerde_groei,  
        model.surface_limited_growth,  
        model.Von_Bertalanffy_groei,  
        model.Gompertz_growth  
    ]  
``  
  
Bij alle modellen moet je de eerste start paramaters instellen. Afhankelijke van de modellen zijn een aantal standaard parameters ingesteld.  
Zie hierboven waar deze parameters aan moeten voldoen.  
  
### Step 2    
Vervolgens zijn er verschillende mogelijkheden.   
  
#### 1. Fitten   
Het fitten van je model, geef een list met x waardes en een list met de passende y waardes. Hierbij worden de beste parameters gezocht    
volgens de direct search volgens Hooke en Jeeves met behulp van de MSE.     
Het is mogelijk om de verschillende solve method mee te geven, standaard word het Runge Kutta gekozen.    
  
`` 
model.fit(x, y, solvingsystem)
``  
  
#### 2. Solve method  
In deze package zitten 2 solved methods, Runge Kutta en Eulers. Geef de begin y waarde mee, de solve berekend de volgend. Hierin kan je een   keuze maken doormiddel van:  
  
``
[  
model.solve_eulers(y0),  
model.solve_rungekutta(y0)  
]  
``
  
### Step 3  
Accuraatheid berekenen van je model doormiddel van BIC.  
`` 
[  
    model.BIC(x, y)
]  
``  
## Getting started  
  
### Requirements:  
- Python 3.x  
- numpy  
- matplotlib  
- pandas  
  
### Clone the repo  
  
``
git clone https://github.com/Larissavf/casus_C  
``
  
This repo contains 2 files:    
    - TGModels.py : De package  
    - product.ipynb : Voorbeeld document met hoe je kan werken met de code  
  
### Run it  
Kies de benodigde stappen die je nodig ben voor jou project.  


