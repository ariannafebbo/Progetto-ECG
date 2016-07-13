
# coding: utf-8

# In[ ]:

import csv
import scipy
import numpy as np
from biosppy.signals import ecg
from termcolor import colored
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

####HO A DISPOSIZIONE 5 ESEMPI DI ECG (ecg1, ecg2, ecg3, ecg4, ecg5)
with open('ecg1.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    
    tempo=list()
    ecgdata=list()
    
    for row in readCSV:
        tempo.append(float(row[0]))
        ecgdata.append(float(row[1]))
        
    T=tempo[2]-tempo[1] #periodo di campionamento
    fs=1/T #frequenza di campionamento


ecg_out = ecg.ecg(ecgdata, fs, show=False)

#distanza temporale tra picchi R, calcolo il tempo in cui si verificano il primo e il secondo picco
t1=tempo[ecg_out['rpeaks'][0]]
t2=tempo[ecg_out['rpeaks'][1]]

#calcolo delta T
deltaT0=(t2-t1)/20 #valore convenzionale
deltaIndice=int(deltaT0/T) #numero di celle di cui mi devo spostare


#calcolo K per picco T
K0 = (t2-t1)*0.6 #valore convenzionale dell'intervallo di tempo
K = int(K0/T) #numero effettivo di celle di cui mi devo spostare 


###ECGDATA FILTRATO
ecgdataF=ecg_out['filtered']


###SCORRO I VARI PICCHI R
indiciQ = list()
indiciT = list()

for i in ecg_out['rpeaks']:
    
    ###CONTROLLO CHE QUANDO ESEGUO IL CICLO FOR NON VADO A LEGGERE CELLE INESISTENTI DEL VETTORE ECGDATA
    #la prima condizione si riferisce a quando il picco R cade all'inizio del vettore, mentre la seconda a quando il picco R cade alla fine
    
    if i-K>0 and i+deltaIndice+K+1<len(ecgdataF):
    #se la condizione è vera si può effettuare la ricerca, altrimenti si passa al picco R successivo
       
        ###RICERCA PICCO T
        massimo = ecgdataF[i+deltaIndice]
        indiceMassimo = i+deltaIndice
        
        for j in range(i+deltaIndice, i+deltaIndice+K+1):
            if ecgdataF[j] > massimo:
                massimo = ecgdataF[j]
                indiceMassimo = j
        indiciT.append(indiceMassimo) #per ogni i viene aggiunto un elemento al vettore che dice dove sono gli indici dei picchi T
        
        ####RICERCA PICCO Q
        minimo = ecgdataF[i-K]
        indiceMinimo = i-K
        
        for j in range(i-2*deltaIndice, i+1):
            if ecgdataF[j] < minimo:
                minimo = ecgdataF[j]
                indiceMinimo = j
        indiciQ.append(indiceMinimo) #per ogni i viene aggiunto un elemento al vettore che dice dove sono gli indici dei picchi Q

#estraggo valori di tempo e ecg in corrispondenza degli indici in cui si trovano i picchi T
tempoT=list()
ecgT=list()

for i in indiciT:
    tempoT.append(tempo[i])
    ecgT.append(ecgdataF[i])

#estraggo valori di tempo e ecg in corrispondenza degli indici in cui si trovano i picchi Q
tempoQ=list()
ecgQ=list()

for i in indiciQ:
    tempoQ.append(tempo[i])
    ecgQ.append(ecgdataF[i])
    
   
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
    
  
######TROVIAMO A PARTIRE DAL PICCO Q QUANDO L'ECG RAGGIUNGE IL LIVELLO DI ZERO, IN MODO DA TROVARE LA DISTANZA INIZIO-PICCO DELL'ONDA Q    

indiciInizioQ=list()

for i in indiciQ:

    j=0
    while ecgdataF[i-j]<0:
        j=j+1

    indiciInizioQ.append(i-j)
            
    
#estraggo valori di tempo e ecg in corrispondenza degli indici in cui si trovano le intersezioni con lo zero prima dei vari picchi Q
tempoInizioQ=list()
ecgInizioQ=list()

for i in indiciInizioQ:
    tempoInizioQ.append(tempo[i])
    ecgInizioQ.append(ecgdataF[i])

    

######TROVIAMO A PARTIRE DAL PICCO T QUANDO L'ECG RAGGIUNGE IL LIVELLO DI ZERO, IN MODO DA TROVARE L'INTERVALLO DI TEMPO PICCO-FINE DELL'ONDA T    

indiciFineT=list()

for i in indiciT:
        
    j=0
    while ecgdataF[i+j]>0:
        j=j+1

    indiciFineT.append(i+j)

            
    
#estraggo valori di tempo e ecg in corrispondenza degli indici in cui si trovano le intersezioni con lo zero dopo i vari picchi T
tempoFineT=list()
ecgFineT=list()

for i in indiciFineT:
    tempoFineT.append(tempo[i])
    ecgFineT.append(ecgdataF[i])
    
    
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

print('Questo programma analizza la lunghezza degli intervalli QT relativi agli ECG di diversi pazienti.\nLa presenza di un solo valore critico viene interpretato come caso patologico.\nSANO: QT < 0.430 s\nBORDERLINE: 0.431 s <= QT <= 0.450 s\nMALATO: QT > 450 s\n')
print('Lunghezza degli intervalli QT:')
QT=list()

#definisco variabile che se il paziente è malato vale 1, altrimenti vale 0
malato=0
#definisco variabile che se il paziente è borderline vale 1, altrimenti vale 0
borderline=0

for i in range(0, len(tempoFineT)):
    valoreQT=tempoFineT[i]-tempoInizioQ[i]
    QT.append(valoreQT)
    
    #BORDERLINE
    if valoreQT>=0.431 and valoreQT<=0.450:
        borderline=1
        print colored(valoreQT, 'yellow')
    
    #MALATO
    elif valoreQT>0.450:
        malato=1
        print colored(valoreQT, 'red')
    
    else:
        print colored(valoreQT, 'green')


if malato==1:
    print('\nPAZIENTE MALATO')
    
elif borderline==1:
    print('\nPAZIENTE BORDERLINE')  

else:
    print('\nPAZIENTE SANO')


####################################################################################    
##########VISUALIZZA I GRAFICI
plt.plot(tempo, ecgdataF, 'b', tempoT, ecgT, 'ro', tempoQ, ecgQ, 'go', tempoInizioQ, ecgInizioQ, 'yo', tempoFineT, ecgFineT, 'ko')
plt.ylabel('ECG Filtrato [mV]')
plt.xlabel('Tempo [s]')

blue_patch = mpatches.Patch(color='blue', label='ECG Filtrato')
red_patch = mpatches.Patch(color='red', label='Picchi T')
green_patch = mpatches.Patch(color='green', label='Picchi Q')
yellow_patch = mpatches.Patch(color='yellow', label='Inizio onda Q')
black_patch = mpatches.Patch(color='black', label='Fine onda T')

plt.legend(handles=[blue_patch, green_patch, red_patch, yellow_patch, black_patch])

plt.grid()
plt.show()



# In[5]:

get_ipython().system(u'pip install termcolor')


# In[ ]:



