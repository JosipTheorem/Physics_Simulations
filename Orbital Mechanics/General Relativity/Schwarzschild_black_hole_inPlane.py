import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import solve_ivp



def Schw_func(lamda,y):
    
    r  = y[0]
    v_r = y[1]
    phi  = y[2]
    v_phi = y[3]  
    
    a_r = (r-3*M) * v_phi**2
    a_phi = -2*v_r*v_phi/r
    
    return [v_r,a_r,v_phi,a_phi]



def stop_at_2M(lamda, y):
    return y[0] - 2*M
stop_at_2M.terminal = True     
stop_at_2M.direction = -1      






if __name__ == '__main__':


    ########## Schwarzschild black hole
    M = 2 
    rays_n = 15

    x_0_values = np.ones(rays_n)*0
    y_0_values = np.linspace(2*M+0.1, 3*M + 2, rays_n)

    vx_values = np.ones(rays_n)*(-1)
    vy_values = np.ones(rays_n)*(0)

    x_0_values = np.concatenate([np.zeros(rays_n), [3*M+0.5, 3*M+0.5]])
    y_0_values = np.concatenate([np.linspace(2*M+0.1, 3*M + 2, rays_n), [0, 0]])
    vx_values  = np.concatenate([np.full(rays_n, -1), [0, 0]])
    vy_values  = np.concatenate([np.zeros(rays_n), [1, -1]])

    lamda = np.linspace(0, 12, 10000)
    total_rays = len(x_0_values)



    lamda_span = (0, 50)

    fig, ax = plt.subplots(figsize=(10, 10))

    for i, (x_0, y_0, vx, vy) in enumerate(zip(x_0_values, y_0_values, vx_values, vy_values)):
        
        r_0 = np.sqrt(x_0**2 + y_0**2)
        phi_0 = np.arctan2(y_0, x_0)

        v_r0 = (x_0 * vx + y_0 * vy) / r_0
        v_phi0 = (x_0 * vy - y_0 * vx) / r_0**2

        ini = [r_0, v_r0, phi_0, v_phi0]

        sol = solve_ivp(Schw_func, lamda_span, ini, events=stop_at_2M, max_step=0.01, rtol=1e-8, atol=1e-8)

        r = sol.y[0]
        phi = sol.y[2]

        x = r * np.cos(phi)
        y = r * np.sin(phi)
        
        color = 'green' if i >= total_rays - 2 else 'blue'
        ax.plot(x, y, color=color, zorder=0)
        #ax.scatter(x[0], y[0], c='black', label='Start')


    x_0 = 3*M-0.5
    y_0 = 0

    vx = 0
    vy = 1

    r_0 = np.sqrt(x_0**2 + y_0**2)
    phi_0 = np.arctan2(y_0, x_0)

    v_r0 = (x_0 * vx + y_0 * vy) / r_0
    v_phi0 = (x_0 * vy - y_0 * vx) / r_0**2

    ini = [r_0, v_r0, phi_0, v_phi0]

    sol = solve_ivp(Schw_func, lamda_span, ini, events=stop_at_2M, max_step=0.01, rtol=1e-8, atol=1e-8)

    r = sol.y[0]
    phi = sol.y[2]

    x = r * np.cos(phi)
    y = r * np.sin(phi)

    ax.plot(x, y, color='blue', zorder=0)
    ax.scatter(3*M+0.5, 0, c='black')





    #Sch. radijus
    filled_circle = plt.Circle((0, 0), 2*M, facecolor='pink', edgecolor='red',zorder=2)
    ax.add_artist(filled_circle)

    #Foton sfera
    theta = np.linspace(0, 2 * np.pi, 150)
    x_circle = 2*M * np.cos(theta)
    y_circle = 2*M * np.sin(theta)
    ax.plot(x_circle, y_circle, color='red', label='r = 2M')


    x_circle = 3*M * np.cos(theta)
    y_circle = 3*M * np.sin(theta)
    ax.plot(x_circle, y_circle, color='orange', label='r = 2M')


    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    #ax.legend()


    plt.show()
