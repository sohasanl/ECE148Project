import sympy as sp
import numpy as np
import itertools
import pandas as pd
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from sympy.core.function import Function

def round_down(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n * multiplier) / multiplier


def GPSPathFinder(initial_location, final_location, obstacles_location_1, obstacles_location_2):
    steps =[0.0, 0.0]
    steps[0] = round_down(((final_location[0] - initial_location[0])/100), 10)
    steps[1] = round_down(((final_location[1] - initial_location[1])/100), 10)

    a1 = [0.0, 0.0]
    a1[0] = round_down(((obstacles_location_1[0] - initial_location[0]) / steps[0]), 10)
    a1[1] = round_down(((obstacles_location_1[1] - initial_location[1]) / steps[1]), 10)
    a2 = [0.0, 0.0]
    a2[0] = round_down(((obstacles_location_2[0] - initial_location[0]) / steps[0]), 10)
    a2[1] = round_down(((obstacles_location_2[1] - initial_location[1]) / steps[1]), 10)

    
    initial_loc = [0.0 , 0.0]
    cur_loc = initial_loc
    
    final_loc = [100.0, 100.0]
    mu = [a1, a2]
    sigma = [[50.0, 0.0], [0.0, 50.0]]
    u, v = sp.symbols('u v')
    w = (np.power(final_loc[0]-u, 2) + np.power(final_loc[1]-v,2))/20000 + \
        1e4 * (1/(sp.det(sp.Matrix(sigma))*2*np.pi) )  *  sp.exp(-(1/2) * np.transpose([(u-mu[0][0]),(v-mu[0][1])]).dot( np.linalg.pinv(sigma).dot([(u-mu[0][0]),(v-mu[0][1])]))) + \
        1e4 * (1/(sp.det(sp.Matrix(sigma))*2*np.pi) )  *  sp.exp(-(1/2) * np.transpose([(u-mu[1][0]),(v-mu[1][1])]).dot( np.linalg.pinv(sigma).dot([(u-mu[1][0]),(v-mu[1][1])])))    

    dx = sp.diff(w, u)
    dy = sp.diff(w, v)
    dz = sp.lambdify([u,v],w,'numpy') 
 

    x = list(np.arange(0.0, 104.0, 4.0))
    y = list(np.arange(0.0, 104.0, 4.0))
    x, y = np.meshgrid(x, y)
    z = dz(x,y)
    

    fig1 = plt.figure()
    ax = fig1.add_subplot(projection='3d')
    ax.plot_wireframe(x, y, z,  rstride=1, cstride=1, alpha = 0.5)
    ax.scatter(initial_loc[0],initial_loc[1], dz(initial_loc[0], initial_loc[1]),c = 'b',marker = '*',alpha = 0.8)

    ax.text(initial_loc[0],initial_loc[1], dz(initial_loc[0], initial_loc[1]), " Source", color='red')
    ax.scatter(final_loc[0],final_loc[1], dz(final_loc[0], final_loc[1]),c = 'r',marker = '*',alpha = 0.8) 
    ax.text(final_loc[0],final_loc[1], dz(final_loc[0], final_loc[1]), " Goal", color='blue')
    
    p = []
    while not(cur_loc == final_loc) :
        if (cur_loc[0] > final_loc[0] or cur_loc[1] > final_loc[1]) :
            break
        dX = round_down(dx.evalf(subs = {u : cur_loc[0], v : cur_loc[1]}), 10)
        dY = round_down(dy.evalf(subs = {u : cur_loc[0], v : cur_loc[1]}), 10)
        d = sp.sqrt(dX**2 + dY**2)
        alpha = round_down(1/d, 10)
        
        cur_loc = list(cur_loc - (alpha * np.array([dX, dY],dtype ='float')))
        cur_loc = [round_down(inx, 10) for inx in cur_loc]
        ax.scatter(cur_loc[0],cur_loc[1],float(w.evalf(subs = {u : cur_loc[0], v : cur_loc[1]})),c = 'g',marker = '*',alpha = 0.8)
        p.append(cur_loc)
        
    for i in p:
        print(i[0])
        i[0] = (i[0] * steps[0]) + initial_location[0]
        i[1] = (i[1] * steps[1]) + initial_location[1]
    plt.show()
    return p
 
 
 
 
 
 
 
    
initial_location = [32.865816, -117.210019]
final_location = [32.865695, -117.210167]  
obstacles_location_1 = [32.865797, -117.210067]
obstacles_location_2 = [32.865731, -117.210106]
path = GPSPathFinder(initial_location, final_location, obstacles_location_1, obstacles_location_2)
path.append(obstacles_location_1)
path.append(obstacles_location_2)
df = pd.DataFrame(path)
df.columns = ['Latitude', 'Longitude']
df.to_csv (r'C:\Users\Omid\Desktop\export_dataframe.csv', index = False, header=True)
#go to GoogleMyMaps and import the .csv file
print(df)


    
