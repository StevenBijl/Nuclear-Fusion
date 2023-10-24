# -*- coding: utf-8 -*-
"""
@author: Steven Bijl
@Project: Nuclear Fusion
"""
import numpy as np
import matplotlib.pyplot as plt
import math





"""
### Seperate the data into energy and MeV/g/cm^2
"""

file = np.loadtxt('Stopping Power & Range.txt')

xaxis_range = file[:, 0]
yaxis_range = file[:, 2]


"""
### Finding the energy deposit in the first and second layer
"""  
### Initial Conditions:

# Traverse distance of the proton in first layer
length_first_layer = 0.03
angle_of_incidence = math.radians(10)
traverse_distance = math.hypot(length_first_layer, length_first_layer *np.sin(angle_of_incidence))
# traverse_distance = 0.03
# print(traverse_distance)

#%%

"""
Here we use the range to approximate the energy deposit of the proton
"""
# desired function (exponential)

polyfit_range = np.polyfit(xaxis_range, yaxis_range, 10)    
trendpoly_range = np.poly1d(polyfit_range)

# curve = trendpoly(xaxis_range)
# print(trendpoly_range(14))


#Verify curve fitting by plotting
fig = plt.figure()
ax = fig.subplots()
ax.scatter(xaxis_range, yaxis_range, color = 'b', s = 5)
ax.set_ylabel('Range [cm]')
ax.set_xlabel('Proton Energy [MeV]')
ax.plot(xaxis_range,trendpoly_range(xaxis_range), color = 'g')
ax.grid()
plt.show()

#%%

"""
Here we use the range to approximate the energy deposit of the proton
"""

A1 = np.array([0,0.1659,0.01049,0.01728])
A2 = np.array([-0.2012,0.2870,0.0378,-0.0011])




def proton_range(energy,mixing_fast,mixing_slow):
    # Curve fitting the Projected Range (Data taken from NIST)
    upper_limit = energy
    lower_limit = energy-5
    
    
    x = np.linspace(lower_limit,upper_limit, 10000)

    polyfit_range = np.polyfit(xaxis_range, yaxis_range, 10)    
    trendpoly_range = np.poly1d(polyfit_range)
    
    range_left = trendpoly_range(energy) - traverse_distance
    fit = trendpoly_range(x)
    idx = (np.abs(fit-range_left)).argmin()

    # print(x[idx])
    
    if x[idx] < 0:
        x[idx] = 0
    
    if x[idx] > 1.85:
        
        F2 = A2[0] + A2[1]*energy + A2[2]* energy**2 + A2[3] * energy**3
        F1 = A2[0] + A2[1]*x[idx] + A2[2]* x[idx]**2 + A2[3] * x[idx]**3
        flayer = F2-F1
        
        # print(1)
        
        # print(energy-x[idx])
        # print(flayer)
        
        # print(x[idx])
        # print(F1)
        
    else:
        if energy > 1.85:
            F2 = A2[0] + A2[1]*energy + A2[2]* energy**2 + A2[3] * energy**3
            F1 = A1[0] + A1[1]*x[idx] + A1[2]* x[idx]**2 + A1[3] * x[idx]**3
            flayer = F2-F1
            
            # print(2)
            
            # print(energy-x[idx])
            # print(flayer)
            
            # print(x[idx])
            # print(F1)
            
        else:
            F2 = A1[0] + A1[1]*energy + A1[2]* energy**2 + A1[3] * energy**3
            F1 = A1[0] + A1[1]*x[idx] + A1[2]* x[idx]**2 + A1[3] * x[idx]**3
            flayer = F2-F1
            
            # print(3)
            
            # print(energy-x[idx])
            # print(flayer)
            
            # print(x[idx])
            # print(F1)
            

    # Energy in the short gate:
    # flayer = energy-x[idx]                          # Initial Proton Energy - Proton Energy after the first layer
    flayer_fast_mixing = flayer*mixing_fast         # Mixing in the first layer
    flayer_fast_eff = flayer_fast_mixing*0.68       # First layer efficiency
    
    flayer_slow_mixing = F1*(1-mixing_slow)     # Mixing in the second layer
    flayer_slow_eff = flayer_slow_mixing*0.41       # Second layer efficiency
    
    short_gate = flayer_fast_eff + flayer_slow_eff  # Energy detected in the short gate
    
    # Energy in the long gate:
    slayer_slow_mixing = F1         # Mixing in the second layer
    slayer_slow_eff = slayer_slow_mixing*0.41       # Second layer efficiency
    
    slayer_fast_mixing = flayer*(1-mixing_fast)     # Mixing in the first layer
    slayer_fast_eff = slayer_fast_mixing*0.68       # First layer efficiency
    
    long_gate = slayer_slow_eff + slayer_fast_eff   # Energy Detected in the long gate

    

    return long_gate, short_gate, x[idx]
    

#%%
"""
Plotting short vs long gate for all energies
"""

def plot_graph_marker_range():
    long_gate =[]
    short_gate = []
    markersx = []
    markersy = []
    for i in range(len(xaxis_range)):
        # print(xaxis_range[i])
        proton_energy = xaxis_range[i]
        a, b,c = proton_range(proton_energy,mixing_fast,mixing_slow)
        if xaxis_range[i] == 3.0 or xaxis_range[i] == 4.0 or xaxis_range[i] == 5.0 or xaxis_range[i] == 6.0 or xaxis_range[i] == 7.0 or xaxis_range[i] == 8.0 or xaxis_range[i] == 9.0 or xaxis_range[i] == 10.0 or xaxis_range[i] == 11.0 or xaxis_range[i] == 12.0 or xaxis_range[i] == 13.0 or xaxis_range[i] == 14.0 or xaxis_range[i] == 15.0 or xaxis_range[i] == 16.0 or xaxis_range[i] == 17.0:
            markersx.append(a)
            markersy.append(b)
            long_gate.append(a)
            short_gate.append(b)
            
        else:
            long_gate.append(a)
            short_gate.append(b)


        text = ["3","4","5","6","7","8","9","10","11","12","13","14","15","16","17"]
        
   
    plt.scatter(markersx,markersy, color = "red")
    for i in range(len(markersx)):    plt.annotate(text[i], (markersx[i]+ 0.1, markersy[i]))
    plt.plot(long_gate, short_gate, color = "red")
    plt.show()


#%%
def scaling_diagonal():
    x = np.linspace(0.00001,0.4,5000)

    slopes = []
    for i in range(len(x)):
        mixing_fast = 1-x[i]
        mixing_slow = 1
        
        dydx1 = proton_range(3,mixing_fast,mixing_slow)
        dydx2 = proton_range(3.1,mixing_fast,mixing_slow)
        
        dx = dydx2[0]-dydx1[0]
        dy = dydx2[1]-dydx1[1]
        
        dydx = dy/dx
        slopes.append(dydx) 
         

        # print(dy/dx)
        
    slopes = np.array(slopes)

    idx = (np.abs(slopes-diagonal_PSD)).argmin()
    # print(idx)
    # print(slopes[idx])
    print(f'The mixing from the fast to the slow layer is: {1-x[idx]}')
    print(f'mixing_fast = {1-x[idx]}')
    return 1- x[idx]


# scaling_diagonal()

#%%

def scaling_angle(mixing_fast):
    x = np.linspace(0.00001,0.6,5000)
    
    
    
    slopes = []
    for i in range(len(x)):
        mixing_slow = 1-x[i]
        
        # svsl = short_gate/long_gate
        # print(short_gate)
        # print(long_gate)
        # print(energy[0], channel4D, svsl)
        
        dx,dy,c = proton_range(14,mixing_fast,mixing_slow)
        dydx = dy/dx
        
        slopes.append(dydx) 
        
        
    channel14D = 137000.0/(485000.0+26777/diagonal_PSD)
    #channel14D = 1.367E+05/(5.043E+05)
    slopes = np.array(slopes)
    
    idx = (np.abs(slopes-channel14D)).argmin()
    # print(idx)
    # print(slopes[idx])
    print(f'The mixing from the slow to fast layer is: {1-x[idx]}')
    print(f'mixing_slow = {1-x[idx]}')
    return 1 - x[idx]

#scaling_angle()


#%%

def scaling_vector_range(island_x,island_y,translation):
    long_gate =[]
    short_gate = []
    non_cont = []
    
    for i in range(len(xaxis_range)):
        proton_energy = xaxis_range[i]
        
        
        a,b,c = proton_range(proton_energy,mixing_fast,mixing_slow)
        
        non_cont.append(c)
        
        long_gate.append(a)
        short_gate.append(b)
    
    
    for i in range(len(non_cont)-1):
        if non_cont[i] == 0.0 and non_cont[i+1] > 0:
            t = i+1
    

    
    shift = translation/diagonal_PSD
    
    long_gate = np.array(long_gate)
    short_gate =np.array(short_gate)
    # print(long_gate[t:])
    
    polyfit = np.polyfit(long_gate[t+10:], short_gate[t+10:], 10)    
    trendpoly = np.poly1d(polyfit)
    
      
    
    x = np.linspace(long_gate[t],long_gate[-1],10000)
    

    d = island_y/(island_x+shift)
    zero_crossing = trendpoly(x) -d*x
    
    
    k = 0
    # print(zero_crossing)
    for i in range(len(zero_crossing)-1):
        if zero_crossing[i] > 0 and zero_crossing[i+1] < 0:
            # print( "crossing between indexes {} and {}".format(i, i+1))
            # print(x[i+1])
            # print(a*np.exp(b*x[i+1])+c)
            # print(d*x[i+1])
            scaling_factor = x[i+1]/(island_x+shift)
            # print(scaling_factor)
            
            # print(trendpoly(x[i]),d*x[i])
            # print(d*(island_x+shift))
            
            k += 1
    
    print(f'Number of zero crossings = {k} \n')
    
    print("translation = {}".format(translation))
    print("island_x = {}".format(island_x))
    print("island_y = {}".format(island_y),'\n')
    print("Your scaling factor is:{}".format(scaling_factor),'\n')     
    
    print(f'plt.hist2d({scaling_factor}*(lg+{shift/diagonal_PSD}),{scaling_factor}*sg, bins=[xbins, ybins], norm=mpl.colors.LogNorm(), cmap=plt.cm.jet)')
    
    
    fig = plt.figure()
    ax = fig.subplots()
    ax.scatter(long_gate, short_gate, color = 'b', s = 5)
    ax.plot(x,trendpoly(x), color = 'g')
    ax.set_ylabel('Short Gate (MeV)')
    ax.set_xlabel('Long Gate (MeV)')
    plt.show()
    plot_graph_marker_range()





#%%

diagonal_PSD = 1.775


### Second approach
proton_energy = 14

mixing_fast = 0.6396069313862772
mixing_slow = 0.8195969273854771

# mixing_fast = scaling_diagonal()
# mixing_slow = scaling_angle(mixing_fast)


# proton_range(proton_energy,mixing_fast,mixing_slow)

# plot_graph_marker_range()


### Getting Coordinates

translation = 26777
island_x=485000.0
island_y=137000.0


scaling_vector_range(island_x,island_y,translation)




