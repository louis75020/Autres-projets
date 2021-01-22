from pylab import *

g=9.80664 #m/s^2

#Schema Euler explicite

#1)
#s1 et s2 seront reutilises pour Verlet

def s1(y,l1,l2,m1,m2):
    a=-m2*l1*(y[3]**2)*sin(y[0]-y[1])*cos(y[0]-y[1])
    a=a+g*m2*sin(y[1])*cos(y[0]-y[1])
    a=a-m2*l2*(y[3]**2)*sin(y[0]-y[1])
    a=a-(m1+m2)*g*sin(y[0])
    a=a/(l1*(m1+m2)-m2*l1*(cos(y[0]-y[1])**2))
    return(a)
    
def s2(y,l1,l2,m1,m2):
    a=m2*l2*(y[3]**2)*sin(y[0]-y[1])*cos(y[0]-y[1])
    a=a+g*(m1+m2)*sin(y[0])*cos(y[0]-y[1])
    a=a+(m1+m2)*l1*(y[3]**2)*sin(y[0]-y[1])
    a=a-g*(m1+m2)*sin(y[1])
    a=a/(l2*(m1+m2)-m2*l2*(cos(y[0]-y[1])**2))
    return(a)
    
def f(y,l1,l2,m1,m2):
    f=zeros(4)
    f[2]=s1(y,l1,l2,m1,m2)
    f[3]=s2(y,l1,l2,m1,m2)
    f[0]=y[2]
    f[1]=y[3]
    return(f)

#2)
def stepEuler(yt,h,l1,l2,m1,m2):
    return(yt+h*f(yt,l1,l2,m1,m2))

#3)
def Euler(y0,P,t,l1,l2,m1,m2): #la boucle est dans une fonction pour eviter de la reexecuter a chaque fois
    i=0
    yt=y0
    T=[]
    while i<=P:
        T.append(yt)
        yt=stepEuler(yt,t/P,l1,l2,m1,m2)#le h du DL a ete defini comme t/P
        i+=1
    return(yt,T)
    
#4)voir fonction test1

#5)    
def plot_angles(T,P,t):#t et P seront reutilises dans tous les graphiques afin de redefinir l echelle de temps cela permet une bonne portabilite
    n=[]#echelle temporelle
    n.append(0)
    for i in range(1,P):
        n.append(i*t/P)
    theta1=[]
    theta2=[]
    for i in range(0,P):
        a=T[i]
        theta1.append(a[0])
        theta2.append(a[1])
    #graphique
    plot(n,theta1,label="pendule 1")
    plot(n,theta2,label="pendule 2")
    xlabel("Temps en s")
    ylabel("angles")
    legend()
    show()
    #plt.savefig("evol-angles.pdf")
    
def evol_vitesse_angulaire(T,P,t):
    n=[]
    n.append(0)
    for i in range(1,P):
        n.append(i*t/P)
    thetapoint1=[]
    thetapoint2=[]
    for i in range(0,P):
        a=T[i]
        thetapoint1.append(a[2])
        thetapoint2.append(a[3])
    #graphique
    plot(n,thetapoint1,label="pendule1")
    plot(n,thetapoint2,label="pendule2")
    xlabel("temps en s")
    ylabel("vitesse angulaire")
    legend()
    show()
    #plt.savefig("evol-vitesse-angulaire.pdf")
    
#b)

def plot_trajectoire(T,l1,l2):#les longueurs des pendules sont reutilises pour obtenir la trajectoire a partir de la position angulaire
    x1=[]
    y1=[]
    x2=[]
    y2=[]
    for i in range(0,1000):
        a=T[i]
        x1.append(l1*sin(a[0]))
        y1.append(-l1*cos(a[0]))
        x2.append(l1*sin(a[0])+l2*sin(a[1]))
        y2.append(-l1*cos(a[0])-l2*cos(a[1]))
    #graphique
    plot(x1,y1,label="pendule1")
    plot(x2,y2,label="pendule2")
    xlabel("x, position longitudinale")
    ylabel("y, position latérale")
    legend()
    show()
    #plt.savefig("trajectoire.pdf")

#6)
#a)       
def energy(y,l1,l2,m1,m2):
    Ec=1/2*(m1+m2)*(y[2]**2)*(l1**2)+1/2*m2*(y[3]**2)*(l2**2)+m2*y[2]*y[3]*l1*l2*cos(y[0]-y[1])
    Ep=-(m1+m2)*g*l1*cos(y[0])-m2*l2*g*cos(y[1])
    return(Ec+Ep)
#b)    
def evol_energy(blabla,y0,l1,l2,m1,m2,P,t):
    T=[]
    if(blabla=="Euler"):a=Euler(y0,P,t,l1,l2,m1,m2)
    if(blabla=="Verlet"):a=Verlet(y0,P,t,l1,l2,m1,m2)
    b=a[1]
    for i in range(0,P):
        y=b[i]
        T.append(energy(y,l1,l2,m1,m2))
    return(T)
#c)
def plot_energy_evol(blabla):
    y0=array([pi/6,pi/6,0,0])#cette fonction sera reutilisee pour Verlet avec les memes conditions initiales
    if(blabla=="Euler"):T=evol_energy("Euler",y0,1,1,1,1,1000,10)
    if(blabla=="Verlet"):T=evol_energy("Verlet",y0,1,1,1,1,1000,10)
    n=[]
    n.append(0)
    for i in range(1,1000):
        n.append(i/100)
    #graphique
    plot(n,T,label="energie au cours du temps")
    xlabel("temps en s")
    ylabel("energie en J")
    legend()
    show()
    #plt.savefig("energie.pdf")

def test1():
    y0=array([pi/2,pi/2,0,0])
    a=Euler(y0,1000,2,1,1,1,1)
    #print(a) #4)
    i=0
    while(i<=1000):
        print(a[1][i])
        i+=100
    print(a[0])
    figure()
    plot_angles(a[1],1000,2)#5a)
    figure()
    evol_vitesse_angulaire(a[1],1000,2)
    figure()
    plot_trajectoire(a[1],1,1)#5b)
#on obtient bien le meme mouvement du pendule que sur la video
    figure()
    plot_energy_evol("Euler") #6c)
#le niveau d energie du systeme est decroissant en valeur absolue au cours du temps: le systeme connait des petes d energie
test1()

#Schema de Verlet
#1)
def stepVerlet(yt,h,l1,l2,m1,m2):
    omega=zeros(4)
    #print(yt[0])
    #print(h)
    y1=yt[0]+h*yt[2]+(h**2)/2*s1(yt,l1,l2,m1,m2)
    y2=yt[1]+h*yt[3]+(h**2)/2*s2(yt,l1,l2,m1,m2)
    omega=array([y1,y2,yt[2],yt[3]])
    y3=yt[2]+h/2*(s1(omega,l1,l2,m1,m2)+s1(yt,l1,l2,m1,m2))
    y4=yt[3]+h/2*(s2(omega,l1,l2,m1,m2)+s2(yt,l1,l2,m1,m2))
    return(array([y1,y2,y3,y4]))
#2)    
def Verlet(y0,P,t,l1,l2,m1,m2):
    T=[]
    yt=y0
    for i in range(0,P):
        T.append(yt)
        yt=stepVerlet(yt,t/P,l1,l2,m1,m2)
    return(yt,T)
    
def test2():
    y0=array([pi/2,pi/2,0,0])
    a=Verlet(y0,1000,2,1,1,1,1)
    figure()
    plot_trajectoire(a[1],1,1)
    figure()
    plot_angles(a[1],1000,2)
    figure()
    evol_vitesse_angulaire(a[1],1000,2)
    figure()
    plot_energy_evol("Verlet")#3)
#malgré de plus fortes variations, le niveau d'energie reste a peu pres constant
test2()

#Chaos
#1)
def plot_espace_phases(T,P): #T tableau de vecteurs représentant les etats successifs de la position aux instants t
    theta1=[]
    theta2=[]
    thetapoint1=[]
    thetapoint2=[]
    for i in range(0,P):
        a=T[i]
        theta1.append(a[0])
        theta2.append(a[1])
        thetapoint1.append(a[2])
        thetapoint2.append(a[3])
    plot(theta1,thetapoint1,label="espaces-phases pendule1")
    plot(theta2,thetapoint2, label="espace-phases pendule2")
    xlabel("angle en rad")
    ylabel("vitesse angulaire en rad/s")
    legend()
    show()
    
def test3():#2a)
    y0=array([10**(-8),sqrt(2)*10**(-8),0,0])
    T=Verlet(y0,1000,20,1,1,1,1)
    print(T[1])
    figure()
    subplot(211)
    plot_trajectoire(T[1],1,1)
    subplot(212)
    plot_espace_phases(T[1],1000)
#la courbe des angles (tethapoint par rapport a theta) est parfaitement sinusoïdale
#les courbes des vitesses angulaires par rapport aux angles sont des ellipses (ici très régulières)
test3()

def test4():#2b)
    y0=array([10**(-8),10**(-8),0,0])
    T=Verlet(y0,1000,20,1,1,1,1)
    figure()
    subplot(211)
    plot_trajectoire(T[1],1,1)
    subplot(212)
    plot_espace_phases(T[1],1000)
#la courbe des angles est moins sinusoidale que pour 2a)
#il en resulte que pour la courbe de la vitesse angulaire par rapport aux angles,il s'agit toujours d'ellipses mais cette fois-ci nettement moins regulieres
test4()

def test5():#3)
    alpha=[10**(-8),5.10**(-3),1,1.2,1.7,1.9]
    for i in range(0,len(alpha)):
        figure()
        y0=array([alpha[i],alpha[i],0,0])
        a=Verlet(y0,1000,20,1,1,1,1)
        subplot(211)
        plot_trajectoire(a[1],1,1)
        subplot(212)
        plot_espace_phases(a[1],1000)
    #pour alpha[0] et alpha[1] les mouvements sont a peu pres sinusoidaux, les vitesses par rapport aux positions sont elliptiques idem que pour le test4, 
    #pour alpha[2] et alpha[3] les mouvements comme la vitesse sont "etonamment" reguliers: on peut aisement definir une periode pour les positions, quant aux vitesses il ne s agit plus d ellipses mais de formes qui se supperposent bien au cours du temps, mais les oscillations sont beaucoup plus fortes que dans le cas precedent
    #pour alpha[4] et alpha[5] les formes des courbes sont assez indescriptibles le mouvement est assez "chaotique" il ne decoule des formes aucune geometrie exploitable, il y a bien une transition vers un mouvement "chaotique" lorsque alpha augmente
test5()