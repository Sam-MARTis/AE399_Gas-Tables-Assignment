# https://sam-martis.github.io/AE399_Gas-Tables-Assignment/

---
Samanth Martis
23B0046
---
Find the code and repository here: [https://github.com/Sam-MARTis/AE399_Gas-Tables-Assignment](https://github.com/Sam-MARTis/AE399_Gas-Tables-Assignment)

**Try the demo**: **[https://sam-martis.github.io/AE399_Gas-Tables-Assignment/](https://sam-martis.github.io/AE399_Gas-Tables-Assignment/)**


Note, the report is in the repository [FinalReport_23B0046.pdf](FinalReport_23B0046.pdf)
![Image of the calculator](image.png)



## Introduction
Various flow quantities are calculated in the calculator, with adjustable inputs for all three classes of calculations - Isentropic flow, normal shock flow, oblique shock flow.
The result is displayed in a single text area to allow for easy copying to clipboard. 

## Code architecture

All logic code is written in Typescript/Javascript with HTML and CSS used to render it. 

For each case, a class is created. Each class has methods that calculate the secondary quantities from the primary quantities. The primary quantities are quantities that are easy to work with and completely describe the system. For `Isentropic` and `Normalshock` classes, the primary quantity is the mach number while in the `ObliqueShock` class they are Mach number and Beta(Wave angle)

Following the main functions of the class to convert the primary quantities, various inverse functions then follow that take in the secondary quantities as input and output the primary quantities. 
Analytical solution is used when possible for the inverses, and iterative solvers(Either Newton-Raphson or simple iteration stepping) when it is not feasible.

The quantities to be displayed are always calculated from the primary quantities. Hence regardless of the input type, it is converted into the primary quantity via inverses and then fed forward to generate the remaining outputs. This may result in small deviations from intended output but they are next to negligible when the tolerance is made stricter.

Finally some supporting functions are created to handle input and output, including setting the global values of precision, tolerance and $\gamma$

All logic code can be found in the `main.ts` file.


## Isentropic Flows

### Stagnation temperature
$$\frac{T_{t}}{T} = 1 + \frac{\gamma + 1}{2} M^2$$
### Stagnation pressure
$$ \frac{P_{t}}{P} = {\left( \frac{T_{t}}{T} \right)}^{\frac{\gamma}{\gamma-1}} $$
### Stagnation density
$$ \frac{\rho_{t}}{\rho} = {\left( \frac{T_{t}}{T} \right)}^{\frac{1}{\gamma-1}} $$
### Stagnation speed of sound
$$ \frac{a_{t}}{a} = {\left( \frac{T_{t}}{T} \right)}^{\frac{1}{2}} $$
### Critical temperature ratio
$$\frac{T^{*}}{T} = \frac{2}{\gamma-1}$$
### Critical pressure ratio
$$\frac{P^{*}}{P} = \left( \frac{T^{*}}{T} \right)^{\frac{\gamma}{\gamma-1}}$$
### Critical density ratio
$$ \frac{\rho^{*}}{\rho} = {\left( \frac{T^{*}}{T} \right)}^{\frac{1}{\gamma-1}} $$
### Critical cross section ratio
$$\frac{A}{A^{*}} = \frac{1}{M}\left\{ \left( \frac{2}{\gamma+1} \right)\left[ 1 + \frac{{\gamma-1}}{2} M^2 \right] \right\}^{\frac{\gamma+1}{2(\gamma-1)}}$$

### Mach wave angle
$$\beta = \sin^{-1}\left( \frac{1}{M_{1}} \right)$$

## Normal Shock Flows

### Downstream mach number
$$M_{2} = \sqrt{ \frac{{(\gamma-1)M_{1}^{2} +2}}{2\gamma M_{1}^{2} - (\gamma-1)} }$$
### Upstream mach number

$$M_{1} = \sqrt{ \frac{{(\gamma-1)M_{2}^2 +2}}{2\gamma M_{2}^{2}- (\gamma-1)}}$$
### Pressure ratio
$$\frac{p_{2}}{p_{1}} = \frac{2\gamma M_{1}^{2} - (\gamma-1)}{\gamma+1}$$

### Density ratio
$$\frac{\rho_{2}}{\rho_{1}} = \frac{(\gamma+1)M_{1}^2}{(\gamma-1)M_{1}^{2}+ 2}$$
### Temperature Ratio
$$\frac{T_{2}}{T_{1}} = \frac{\frac{P_{2}}{P_{1}}}{\frac{\rho_{2}}{\rho_{1}}}$$

### Sound speed ratio
$$\frac{a_{2}}{a_{1}} = \sqrt{ \frac{T_{2}}{T_{1}} }$$
### Total pressure ratio
$$\frac{P_{2t}}{P_{1t}} = \left\{ \frac{(\gamma+1)}{2} \frac{M_{1}^{2}}{\left( 1+ \frac{{\gamma-1}}{2}M_{1}^{2}  \right)}\right\}^{\frac{\gamma}{\gamma-1}} \left\{ \left( \frac{2\gamma}{\gamma+1} \right)M_{1}^{2}-\left( \frac{{\gamma-1}}{\gamma+1} \right) \right\}^{- \frac{1}{\gamma-1}}$$

### Upstream static vs downstream total
$$\frac{P_{1}}{P_{2t}} = \frac{\frac{P_{1}}{P_{2}}}{\frac{P_{2t}}{P_{2}}}$$

### Mach number from pressure ratio
$$M_1 \;=\; \sqrt{\tfrac{(\gamma - 1) + \tfrac{P_2}{P_1}(\gamma + 1)}{2\gamma}}
$$

### Mach number from density ratio
$$M_1 \;=\; \sqrt{\tfrac{2 \tfrac{\rho_2}{\rho_1}}{(\gamma+1) - \tfrac{\rho_2}{\rho_1}(\gamma - 1)}}
$$
### Mach number from temperature ratio
$$M_1 \;=\; \sqrt{\tfrac{-b + \sqrt{b^2 - 4ac}}{2a}}, \quad a = 2\gamma(\gamma-1),\; b = 4\gamma - (\gamma-1)^2 - (\gamma+1)^2 \tfrac{T_2}{T_1},\; c = -2(\gamma - 1)
$$




## Oblique Shocks

### Deflection angle  
$$
\delta \;=\; \arctan \!\left( \frac{ 2 \,\tfrac{\cos\beta}{\sin\beta} \,\big( (M_1 \sin\beta)^2 - 1 \big) }{ 2 + M_1^2 \big( \gamma + \cos(2\beta) \big) } \right)
$$  

### Downstream Mach number  
$$
M_2 \;=\; \sqrt{ \frac{ (M_1^2 \sin^2 \beta) + \tfrac{2}{\gamma - 1} }{ \left( \tfrac{2\gamma}{\gamma - 1}(M_1^2 \sin^2 \beta) - 1 \right) \sin^2(\beta - \delta) } }
$$  

### Pressure ratio  
$$
\frac{P_2}{P_1} \;=\; \frac{2\gamma (M_1^2 \sin^2 \beta) - (\gamma - 1)}{\gamma + 1}
$$  

### Density ratio  
$$
\frac{\rho_2}{\rho_1} \;=\; \frac{(\gamma + 1)(M_1^2 \sin^2 \beta)}{(\gamma - 1)(M_1^2 \sin^2 \beta) + 2}
$$  
### Temperature ratio  
$$
\frac{T_2}{T_1} \;=\; \frac{\tfrac{P_2}{P_1}}{\tfrac{\rho_2}{\rho_1}}
$$  

### Maximum deflection angle  
$$
\beta_{\delta\max} \;=\; \arcsin \!\left( \sqrt{ \tfrac{\gamma+1}{4\gamma} \;-\; \tfrac{1}{\gamma M_1^2} \Big(1 - \sqrt{ (\gamma+1)\Big(1 + \tfrac{(\gamma-1)}{2} M_1^2 + \tfrac{(\gamma+1)}{16} M_1^4 \Big)} \Big)} \right)
$$  $$
\delta_{max} \;=\; \arctan \!\left( \frac{ 2 \,\tfrac{\cos\beta_{max}}{\sin\beta_{max}} \,\big( (M_1 \sin\beta_{max})^{2} - 1 \big) }{ 2 + M_1^2 \big( \gamma + \cos(2\beta_{max}) \big) } \right)
$$



### Strong and Weak Oblique Shock Solutions

#### Strong & Weak Shock Solutions  
For an upstream Mach number $M_1$, ratio of specific heats $\gamma$, and a deflection angle $\delta$:  

- Let $\beta_w$ = weak shock angle solution  
- Let $\beta_s$ = strong shock angle solution  

Then the solutions are:  

#### Shock angles:
$$
\beta = [\,\beta_s,\; \beta_w\,]
$$  

#### Downstream Mach numbers:
$$
M_2 = \Big[\, M_{2,s} \;=\; M_2(M_1, \beta_s, \gamma),\;\; M_{2,w} \;=\; M_2(M_1, \beta_w, \gamma) \,\Big]
$$  

#### Pressure ratios:  
$$
\frac{P_2}{P_1} = \Bigg[\, \frac{P_2}{P_1}(M_1,\beta_s,\gamma),\;\; \frac{P_2}{P_1}(M_1,\beta_w,\gamma) \,\Bigg]
$$  

#### Density ratios:  
$$
\frac{\rho_2}{\rho_1} = \Bigg[\, \frac{\rho_2}{\rho_1}(M_1,\beta_s,\gamma),\;\; \frac{\rho_2}{\rho_1}(M_1,\beta_w,\gamma) \,\Bigg]
$$  

#### Temperature ratios:  
$$
\frac{T_2}{T_1} = \Bigg[\, \frac{T_2}{T_1}(M_1,\beta_s,\gamma),\;\; \frac{T_2}{T_1}(M_1,\beta_w,\gamma) \,\Bigg]
$$  

---

#### Normal Mach numbers  

The normal Mach numbers corresponding to the weak and strong solutions are:  

##### Upstream normal Machs:  
$$
M_{1n} = \Big[\, M_1 \sin\beta_s,\;\; M_1 \sin\beta_w \,\Big]
$$  

##### Downstream normal Machs:  
$$
M_{2n} = \Big[\, M_2(M_1,\beta_s,\gamma)\,\sin(\beta_s - \delta),\;\; M_2(M_1,\beta_w,\gamma)\,\sin(\beta_w - \delta) \,\Big] 
$$  
