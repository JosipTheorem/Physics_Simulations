import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation 
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.image as mpimg
from scipy.integrate import odeint
from IPython.display import HTML
from matplotlib import cm



import matplotlib
matplotlib.use('TkAgg')
from matplotlib.widgets import Button, Slider


def Potential(r):
    return -k/r

def Fx(x,y):
    r = np.sqrt(x**2 + y**2)
    F_x = -(k*x)/(r**3)
    return F_x

def Fy(x,y):
    r = np.sqrt(x**2 + y**2)
    F_y = -(k*y)/(r**3)
    return F_y

def Energija(x,vx,y,vy):
    r = np.sqrt(x**2 + y**2)
    T = (vx**2 + vy**2)*m/2
    V = -k/r
    return V + T
    
#### Three Point formula for Potential (force)
#def Potential_n(x,y):
#    r = np.sqrt(x**2 + y**2)
#    return -k/r
    
#h = 0.00001
#def Fx_n(x,y):
#    x = 1; y = 2
#    return -(Potential_n(x+h,y)-Potential_n(x-h,y))/(2*h)


def update(val):
    
    #solve again each time 
    T_max = Time.val
    dt = 0.0001
    t = np.arange(0, T_max+dt, dt)
    
    m = 1.8 # Mass of an object
    k = K.val # Potential constant

    r0 = 2.5 # Initial radius
    v0 = V_0.val # Initial velocity
    alpha0 = Alpha.val # Initial angle of velocity in respect to x-axis
    alpha = alpha0*np.pi/180

    # Coordinates of velocity
    v0x = v0*np.cos(alpha)
    v0y = v0*np.sin(alpha)
    # Initial Cartesian Coordinates
    x0 = X_0.val
    y0 = Y_0.val

    P0 = [x0,v0x,y0,v0y]

    
    solution = odeint(func,P0,t)
    x = solution[:,0]
    y = solution[:,2]
    
    l.set_data(x, y)
    plt.draw()

def func(P,t):
        x=P[0]
        vx = P[1]
        y=P[2]
        vy = P[3]
        return[vx,Fx(x,y)/m,vy,Fy(x,y)/m]


def animate(i):

    current_frame = i * frame_step
    #current_frame = i
    


    time_text.set_text(f'Time: {t[current_frame]:.2f} s')
    energy = Energija(solution[current_frame, 0], solution[current_frame, 1], solution[current_frame, 2], solution[current_frame, 3])
    energy_text.set_text(f'Energy: {energy:.2f} J')
    


    color = colormap(current_frame / len(t))  # Color based on the current frame
    l.set_data(solution[:current_frame, 0], solution[:current_frame, 2])
    l.set_color(color)  # Update trajectory line color
    point.set_data([solution[current_frame, 0]], [solution[current_frame, 2]])  # Update moving point


def update(val):
        
    #solve again each time 
    T_max = Time.val
    dt = 0.0001
    t = np.arange(0, T_max+dt, dt)
    
    m = 1.8 # Mass of an object
    k = K.val # Potential constant

    r0 = 2.5 # Initial radius
    v0 = V_0.val # Initial velocity
    alpha0 = Alpha.val # Initial angle of velocity in respect to x-axis
    alpha = alpha0*np.pi/180

    # Coordinates of velocity
    v0x = v0*np.cos(alpha)
    v0y = v0*np.sin(alpha)
    # Initial Cartesian Coordinates
    x0 = X_0.val
    y0 = Y_0.val

    P0 = [x0,v0x,y0,v0y]

    
    solution = odeint(func,P0,t)
    x = solution[:,0]
    y = solution[:,2]
    
    l.set_data(x, y)
    plt.draw()
    
    



if __name__ == '__main__':


    m = 1.8 # Mass of an object
    k = 1.5 # Potential constant

    r0 = 2.5 # Initial radius
    v0 = 0.5 # Initial velocity
    alpha0 = 45 # Initial angle of velocity in respect to x-axis
    alpha = alpha0*np.pi/180

    # Coordinates of velocity
    v0x = v0*np.cos(alpha)
    v0y = v0*np.sin(alpha)
    # Initial Cartesian Coordinates
    x0 = r0
    y0 = 0.0
    # Potential
    V = -k/r0
    # Energy
    E = m*v0**2/2 - k/r0



    #Solving
    T_max = 20
    dt = 0.0001
    t = np.arange(0, T_max+dt, dt)
    P0 = [x0,v0x,y0,v0y] #Initial conditions

    ##Solver
    solution = odeint(func,P0,t)

    x = solution[:,0]
    y = solution[:,2]
    vx = solution[:,1]
    vy = solution[:,3]

 ############## Energy
    E_t = Energija(solution[:,0],solution[:,1],solution[:,2],solution[:,3])
    plt.plot(t, E_t)
    plt.ylim(E+0.1,E-0.1)
    plt.xlabel('t')

    plt.ylabel('E')
    plt.show()
 ############## Energy

 ############## video


    plt.rcParams.update({'font.size': 12, 'font.family': 'sans-serif'})  # Use a clean sans-serif font

    fig, ax = plt.subplots()
    ax.set_facecolor('black')
    fig.patch.set_facecolor('black')
    ax.tick_params(axis='both', colors='white')  # Axis ticks and labels
    ax.xaxis.label.set_color('white')  # X-axis label color
    ax.yaxis.label.set_color('white')  # Y-axis label color
    #ax.grid(True, color='white')  # Grid color
    ax.axis([-2,4,-1,4]) # Axis boundaries

    l, = ax.plot([],[], color='white', alpha=0.75) 
    ax.scatter(0,0)

    max_frames = round(len(t)/1000)  # Number of frames you want to show
    frame_step =  1000

    point, = ax.plot([], [], 'wo')



    time_text = ax.text(0.02, 0.95, f'Time: {t[0]:.2f} s', transform=ax.transAxes, color='white', fontsize=12)
    energy_text = ax.text(0.02, 0.90, f'Energy: {E:.2f} J', transform=ax.transAxes, color='white', fontsize=12)



    colormap = cm.plasma

    ax.set_title('Object in Central Force Field', fontsize=16, color='white', weight='bold')
    ax.set_xlabel('X Position', fontsize=12, color='white')
    ax.set_ylabel('Y Position', fontsize=12, color='white')

    # Load the Sun image and add it at the origin (0, 0)
    #sun_img = mpimg.imread('sun.jpg')  # Replace 'sun_image.png' with the path to your Sun image
    #imagebox = OffsetImage(sun_img, zoom=0.02)  # Adjust `zoom` for the image size
    #ab = AnnotationBbox(imagebox, (0, 0), frameon=False)  # Place at the origin (0, 0)
    #ax.add_artist(ab)


        

    ani = matplotlib.animation.FuncAnimation(fig, animate, frames=max_frames, interval=10)
    plt.show()
 ############## video

 ############## slider


    fig = plt.figure(figsize=(10, 6))
    fig.subplots_adjust(left=0.25, bottom=0.4)
    ax = plt.subplot(111)

    l, = plt.plot(x, y, lw=2, color='red')
    plt.axis([-5, 5, -5, 5])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid()

    axn1 = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor='#e4e4e4')
    axn2 = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor='#e4e4e4')
    axn3 = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor='#e4e4e4')
    axn4 = plt.axes([0.25, 0.2, 0.65, 0.03], facecolor='#e4e4e4')
    axn5 = plt.axes([0.25, 0.25, 0.65, 0.03], facecolor='#e4e4e4')
    axn6 = plt.axes([0.25, 0.3, 0.65, 0.03], facecolor='#e4e4e4')


    Alpha = Slider(axn1, 'alpha_0', 0, 180, valinit=20)
    V_0 = Slider(axn2, 'v_0', 0, 5, valinit=0.5)
    X_0 = Slider(axn3, 'x_0', -10, 10, valinit=2.5)
    Y_0 = Slider(axn4, 'y_0', -10, 10, valinit=0)
    K = Slider(axn5, 'k', 0.1, 10, valinit=1.5)
    Time = Slider(axn6, 'Time', 1, 100, valinit=20)


        
    Alpha.on_changed(update)
    V_0.on_changed(update)
    X_0.on_changed(update)
    Y_0.on_changed(update)
    K.on_changed(update)
    Time.on_changed(update)


    plt.show(block=True)

 ############## slider

