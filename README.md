### User's guide

#### Step 1  
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

#### Step 2  
Vervolgens zijn er verschillende mogelijkheden. 

##### 1. Fitten 
Het fitten van je model, geef een list met x waardes en een list met de passende y waardes. Hierbij worden de beste parameters gezocht  
volgens de direct search volgens Hooke en Jeeves met behulp van de MSE.   
Het is mogelijk om de verschillende solve method mee te geven, standaard kiest het Runge Kutta.
`` 
model.fit(x, y, solvingsystem)
``  

##### 2. Solve method
In deze package zitten 2 solved methods, Runge Kutta en Eulers. Geef de begin y waarde mee, de solve berekend de volgend. Hierin kan je een keuze maken doormiddel van:
``
[
model.solve_eulers(y0),
model.solve_rungekutta(y0)
]
``
#### Step 3
Accuraatheid berekenen van je model doormiddel van BIC.
`` 
[
    model.BIC(x, y)
]
``
### Getting started

#### Requirements:
- Python 3.x
- numpy
- matplotlib
- pandas

#### Clone the repo
``
git clone https://github.com/Larissavf/casus_C
``

This repo contains 2 files:  
    - TGModels.py : De package
    - product.ipynb : Voorbeeld document met hoe je kan werken met de code

#### Run it
Kies de benodigde stappen die je nodig ben voor jou project.


