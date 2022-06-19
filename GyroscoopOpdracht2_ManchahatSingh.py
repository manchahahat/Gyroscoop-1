#!/usr/bin/env python
# coding: utf-8

# In[15]:


#importeren van modules
import numpy as np
import matplotlib.pyplot as plt
import gissa2 as gs


# In[16]:


#stapgrootte en aantal stappen
h = 20000
dt = 0.000005

#parameters stap 1
v0 = 15
vdrive = 1.5
e0 = 8.85e-12
l_drive = 0.2e-3
d_drive = 2e-6
w = 3e-6
N_drive = 100

#parameters stap 2
m_drive = 4.2818e-9
k_drive = 1.0565
gamma_drive = 1.9217e-6
pi = np.pi

#parameters stap 3
ohm = 10/60

#parameters stap 4
m_sense = 5.1328-9
k_sense = 5.0712
gamma_sense = 3.2284e-06

#parameters stap 5
l_sense = 0.2e-3
d_sense = 2e-6
N_sense = 40
Vdc = 15

#lijst aanmaker
def listmaker(h):
    y = np.zeros(h)
    y[0] = 0
    return y


# In[17]:


#resonantiefrequentie defineren
f = 1/(2*pi*(m_drive/k_drive)**0.5)
f


# In[18]:


#tijd defineren
t = listmaker(h)
for i in range(h-1):
    t[i+1] = t[i] + dt


# In[19]:


#spanning over de condensatoren in de drive mode definiëren
vt= listmaker(h)
for i in range (h):
    vt[i] = v0 + vdrive *(np.cos(2 * pi * f * t[i]))


# In[20]:


#capaciteit definitie
C = 8.85E-12 *((l_drive*w)/d_drive)


# In[21]:


#elektrische kracht definitie
Fel = N_drive * 0.5 *vt**2 *C


# In[22]:


#constantes van de differentiaal vergelijking voor het massaveersysteem in de drive- en sensemode
a_xdrive = (k_drive-2*m_drive/(dt**2))/(m_drive/(dt**2)+gamma_drive/(2*dt))
b_xdrive = (m_drive/(dt**2)-gamma_drive/(2*dt))/(m_drive/(dt**2)+gamma_drive/(2*dt))
c_xdrive = Fel/(m_drive/(dt**2)+gamma_drive/(2*dt))

a_xsense = dt**2/(m_sense+dt*gamma_sense)
b_xsense = (2*m_sense+dt*gamma_sense-k_sense*dt**2)/(m_sense+dt*gamma_sense)
c_xsense = m_sense/(m_sense+dt*gamma_sense)


# In[23]:


#uitwijking in de drive mode definitie
x_drive = listmaker(h)

for i in range(h-1):
    x_drive[i+1] = -a_xdrive*x_drive[i] - b_xdrive*x_drive[i-1] + c_xdrive[i]*np.sin(2*pi*f*t[i+1])


#snelheid in de drive mode definitie
v_xdrive = listmaker(h)
for i in range(h-1):
    v_xdrive[i+1] = (x_drive[i+1]-x_drive[i])/dt


# In[24]:


#absolute waarde nemen van de snelheid in de drive mode
vabs = abs(v_xdrive)

#definitie van de coriolis kracht
x_sense = listmaker(h)
Fcori = listmaker(h)
for i in range(h):
    Fcori[i] = -2*m_drive*ohm*v_xdrive[i]
    
#uitwijking in de sense mode definitie
for i in range(h-2):
    x_sense[i+2] = a_xsense*Fcori[i]+b_xsense*x_sense[i+1]-c_xsense*x_sense[i]


# In[25]:


#definitie van de lading, stroomsterkte, capaciteit en de aflgeleide van de capaciteit in de sense mode
Q = listmaker(h)
I = listmaker(h)
C = listmaker(h)
C2 = listmaker(h)
for i in range(h):
    C[i] = N_sense*e0*((l_sense*w)/(d_sense+x_sense[i]))

for i in range(h-1):
    C2[i+1]= (C[i+1]-C[i])/dt
    
for i in range(h):
    Q[i] = C2[i]*Vdc
    
for i in range(h-1):
    I[i+1] = (Q[i+1]-Q[i])/dt


# In[26]:


#het plotten van de grafieken
plot1 = gs.plot(t, vt,'t', 'V', 's', 'V', 0, 0.1)
plot2 = gs.plot(t, Fel, 't', 'F', 's', 'N', 0, 0.1)
plot3 = gs.plot(t, x_drive, 't', 'x', 's', 'm', 0, 0.1)
plot4 = gs.plot(t, Fcori, 't', 'F', 's', 'N', 0, 0.1)
plot5 = gs.plot(t, x_sense, 't', 'x', 's', 'm', 0, 0.1)
plot6 = gs.plot(t, I, 't', 'I', 's', 'A', 0, 0.1)


# In[27]:


#vergrotingsfactor per stap definitie
n1 = max(Fel)/max(vt)
n2 = max(x_drive)/max(Fel)
n3 = max(Fcori)/max(x_drive)
n4 = max(x_sense)/max(Fcori)
n5 = max(I)/max(x_sense)
print(n1,n2,n3,n4,n5)


# In[28]:


#vinden en printen van de in- en outputsignalen van de gyroscoop en de vergelijking met andere referentiewaardes
v_max = max(vt)
Fel_max = max(Fel)
x_drive_max= max(x_drive)
Fcori_max= max(Fcori)
x_sense_max= max(x_sense)
I_max= max(I)

print('v_max       =', v_max, 'V','                       (dit is ongeveer 15 keer kleiner dan wat er thuis uit je stopcontact komt)')
print('Fel_max     =', Fel_max, 'N', '      (dit is ongeveer het zelfde als de zwartekracht die op 4 gemiddelde menselijke cellen wordt uitgeoefend )')
print('x_drive_max =', x_drive_max, 'm', '     (dit is ongeveer 3 keer kleiner dan de diameter van een ribosoom in een menselijke cel)')
print('Fcori_max   =', Fcori_max, 'N', '     (dit is ongeveer 6 keer groter dan de zwaartekraht die op een E-coli bacterie wordt uitgeoefend)')
print('x_sense_max =', x_sense_max, 'm', '      (dit is ongeveer 10 miljard keer kleiner dan de diameter van één haar)')
print('I_max       =', I_max, 'A', '      (dit is ongeveer 6,6 miljard keer kleiner dan de dodelijke hoeveelheid voor een mens)')


# In[ ]:





# In[ ]:




