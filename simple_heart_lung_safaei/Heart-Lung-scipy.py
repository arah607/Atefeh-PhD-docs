import math
from math import cos
from math import floor

# Import numpy
from typing import Optional, Tuple, Iterable
import numpy as np
import matplotlib.pyplot as plt
from numpy import ndarray
from switch import Switch
from math import pi
from math import exp
from math import sqrt
from scipy.integrate import odeint

class Heart_Lung:
    def __init__(self):

        self.T_vc = 0.34  # The duration of ventricles contraction
        self.T_vr = 0.15  # The duration of ventricles relaxation

        self.t_ar = 0.97  # The time when the atria start to relax
        self.T_ar = 0.17  # The duration of atria relaxation
        self.t_ac = 0.80  # The time when the atria start to contraction
        self.T_ac = 0.17  # The duration of atria contraction

        # blood_pressure_Atria_ventricles
        self.E_ra_A = 7.998e+6  # amplitude value of the RA elastance
        self.E_ra_B = 9.331e+6  # baseline value of the RA elastance
        self.V_ra = 20.0e-6  # Initial blood volume of RA
        self.V_ra_0 = 4.0e-6  # dead blood volume of RA

        self.E_rv_A = 73.315e+6  # amplitude value of the RV elastance
        self.E_rv_B = 6.665e+6  # baseline value of the RV elastance
        self.V_rv = 500.0e-6  # Initial blood volume of RV
        self.V_rv_0 = 10.0e-6  # dead blood volume of RV

        self.E_la_A = 9.331e+6  # amplitude value of the LA elastance
        self.E_la_B = 11.997e+6  # baseline value of the LA elastance
        self.V_la = 20.0e-6  # blood volume of LA
        self.V_la_0 = 4.0e-6  # dead blood volume of LA

        self.E_lv_A = 366.575e+6  # amplitude value of the LV elastance
        self.E_lv_B = 10.664e+6  # baseline value of the LV elastance
        self.V_lv = 500.0e-6  # blood volume of LV
        self.V_lv_0 = 5.0e-6  # dead blood volume of LV

        # blood_flow_atria_ventricles
        self.CQ_trv = 34.6427e-6  # triscupid valve coefficient
        self.CQ_puv = 30.3124e-6  # pulmonary valve coefficient
        self.P_pulmonary_artery = 4000.0  # pulmonary arteries
        self.CQ_miv = 34.6427e-6  # mitral valve coefficient
        self.CQ_aov = 30.3124e-6  # aortic valve coefficient
        self.P_root = 0 #blood pressure in the aortic root


        # blood_volume_Atria_ventricles
        self.Q_sup_venacava = 4.21324067213557e-5  # blood flow superior vena cava
        # self.Q_sup_venacava = 0  # blood flow superior vena cava
        self.Q_inf_venacava = 4.11763993353109e-5  # blood flow inferior vena cava
        # self.Q_inf_venacava = 0  # blood flow inferior vena cava
        self.Q_pulmonary_vein = 0  # vein flow

        # pulmonary circulation
        self.Q_pulmonary_artery = 0  # artery flow
        self.C_pulmonary_artery = 0.0309077e-6  # artery compliance
        self.I_pulmonary_artery = 1.0e-6  # artery inductance
        self.R_pulmonary_artery = 10.664e+6  # artery resistance
        self.C_pulmonary_vein = 0.60015e-6  # vein compliance
        self.P_pulmonary_vein = 0  # # vein pressure
        self.R_pulmonary_vein = 1.333e+6  # vein resistance
        self.I_pulmonary_vein = 1.0e-6  # vein inductance

        self.T = 1  # duration of a cardiac cycle


    def define_t_array(self, a, b, num_int):
        ''' This function considers the time as an array '''
        self.t = np.linspace(a, b, num_int)

    def activation_ventricles(self, mt):
        '''  This function calculates the activation ventricles '''
        if (mt >= 0) and (mt <= self.T_vc * self.T):
            self.e_v = 0.5 * (1 - cos(pi * mt / (self.T_vc * self.T)))
        elif (mt > self.T_vc * self.T) and (mt <= (self.T_vc + self.T_vr) * self.T):
            self.e_v = 0.5 * (1 + cos(pi * (mt - self.T_vc * self.T) / (self.T_vr * self.T)))
        elif (mt > (self.T_vc + self.T_vr) * self.T) and (mt < self.T):
            self.e_v = 0

        return self.e_v

    def activation_atria(self, mt):
        ''' This function calculates the  activation atria '''
        if (mt >= 0) and (mt <= (self.t_ar + self.T_ar) * self.T - self.T):
            self.e_a = 0.5 * (1 + cos(pi * (mt + self.T - self.t_ar * self.T) / (self.T_ar * self.T)))
        elif (mt > (self.t_ar + self.T_ar) * self.T - self.T) and (mt <= self.t_ac * self.T):
            self.e_a = 0
        elif (mt > self.t_ac * self.T) and (mt <= (self.t_ac + self.T_ac) * self.T):
            self.e_a = 0.5 * (1 - cos(pi * (mt - self.t_ac * self.T) / (self.T_ac * self.T)))
        elif (mt > (self.t_ac + self.T_ac) * self.T) and (mt <= self.T):
            self.e_a = 0.5 * (1 + cos(pi * (mt - self.t_ar * self.T) / (self.T_ar * self.T)))

        return self.e_a

    def Blood_Press_Atria_Ventricles(self, V_ra, V_rv, V_la, V_lv):
        ''' This function calculates the blood pressure in atria and ventricles '''

        self.P_ra = (self.e_a * self.E_ra_A + self.E_ra_B) * (V_ra - self.V_ra_0)

        # RV pressure
        self.P_rv = (self.e_v * self.E_rv_A + self.E_rv_B) * (V_rv - self.V_rv_0)

        # LA pressure
        self.P_la = (self.e_a * self.E_la_A + self.E_la_B) * (V_la - self.V_la_0)

        # LV pressure
        self.P_lv = (self.e_v * self.E_lv_A + self.E_lv_B) * (V_lv - self.V_lv_0)

        return (self.P_ra, self.P_rv, self.P_la, self.P_lv)

    def Blood_Flow_Atria_ventricles(self):
        '''   This function calculates the blood flow in atria and ventricles '''

        if self.P_ra >= self.P_rv:
            self.Q_ra = self.CQ_trv * sqrt(self.P_ra - self.P_rv)

        else:
            self.Q_ra = 0

        # RV blood flow
        if self.P_rv >= self.P_pulmonary_artery:  # self.P_pulmonary_artery
            self.Q_rv = self.CQ_puv * sqrt(self.P_rv - self.P_pulmonary_artery)
        else:
            self.Q_rv = 0

        # LA blood flow
        if self.P_la >= self.P_lv:
            self.Q_la = self.CQ_miv * sqrt(self.P_la - self.P_lv)
        else:
            self.Q_la = 0

        # LV blood flow
        if self.P_lv >= self.P_root:
            self.Q_lv = self.CQ_aov * sqrt(self.P_lv - self.P_root)
        else:
            self.Q_lv = 0
        return (self.Q_ra, self.Q_rv, self.Q_la, self.Q_lv)


    def ODE_blood_volume_Right_Atria_Ventricles(self, y, t):

        # Store the values in y array
        V_ra = y[0]
        V_rv = y[1]

        # Define the ode formulation for right atrium
        if self.P_ra >= self.P_rv:
            self.der_volume_RA = self.Q_sup_venacava + self.Q_inf_venacava - (self.CQ_trv * sqrt(abs(
                    ((self.e_a * self.E_ra_A + self.E_ra_B) * (V_ra - self.V_ra_0)) - (
                    (self.e_v * self.E_rv_A + self.E_rv_B) * (V_rv - self.V_rv_0)))))


        else:
            self.der_volume_RA = self.Q_sup_venacava + self.Q_inf_venacava





        # Define ode formulation for right ventricle
        if self.P_rv >= self.P_pulmonary_artery:

            self.der_volume_RV = (self.CQ_trv * sqrt(abs(((self.e_a * self.E_ra_A + self.E_ra_B) * (
                        V_ra - self.V_ra_0)) - ((self.e_v * self.E_rv_A + self.E_rv_B) * (V_rv - self.V_rv_0))))) - (
                                         self.CQ_puv * sqrt(abs(((self.e_v * self.E_rv_A + self.E_rv_B) * (
                                         V_rv - self.V_rv_0)) - self.P_pulmonary_artery)))

        else:
            self.der_volume_RV = self.CQ_trv * sqrt(
                abs(((self.e_a * self.E_ra_A + self.E_ra_B) * (V_ra - self.V_ra_0)) - (
                        (self.e_v * self.E_rv_A + self.E_rv_B) * (V_rv - self.V_rv_0))))



        return (self.der_volume_RA, self.der_volume_RV)


    def Solve_ODE_blood_volume_Right_Atria_Ventricles(self, a, b):
        # Initial values for v_ra and v_rv
        y0 = [self.V_ra, self.V_rv]
        # Time between t[i-1] and t[i]
        t = [a, b]
        # Solve ode
        y = odeint(self.ODE_blood_volume_Right_Atria_Ventricles, y0, t)
        # Store values after solve ode
        self.V_ra = y[1,0]
        self.V_rv = y[1,1]

        return (self.V_ra, self.V_rv)

    def ODE_blood_volume_Left_Atria_Ventricles(self, z, t):

        # Store the values in z array
        V_la = z[0]
        V_lv = z[1]

        # Define the ode formulation for left atrium
        if self.P_la >= self.P_lv:
            self.der_volume_LA = self.Q_pulmonary_vein - (self.CQ_miv * sqrt(
                abs(((self.e_a * self.E_la_A + self.E_la_B) * (V_la - self.V_la_0)) - (
                            (self.e_v * self.E_lv_A + self.E_lv_B) * (V_lv - self.V_lv_0)))))

        else:
            self.der_volume_LA = self.Q_pulmonary_vein




        # Define ode formulation for left ventricle
        if self.P_lv >= self.P_root:
            # if self.P_lv >= self.P_root:
            self.der_volume_LV = (self.CQ_miv * sqrt(abs(((self.e_a * self.E_la_A + self.E_la_B) * (
                        V_la - self.V_la_0)) - ((self.e_v * self.E_lv_A + self.E_lv_B) * (V_lv - self.V_lv_0))))) - (
                                             self.CQ_aov * sqrt(abs(((self.e_v * self.E_lv_A + self.E_lv_B) * (
                                                 V_lv - self.V_lv_0)) - self.P_root)))


        else:
            self.der_volume_LV = (self.CQ_miv * sqrt(abs(((self.e_a * self.E_la_A + self.E_la_B) * (
                        V_la - self.V_la_0)) - ((self.e_v * self.E_lv_A + self.E_lv_B) * (V_lv - self.V_lv_0)))))

        # der_volume[0] = self.der_volume_LA
        # der_volume[1] = self.der_volume_LV

        return (self.der_volume_LA, self.der_volume_LV)

    def Solve_ODE_blood_volume_Left_Atria_Ventricles(self, a, b):
        # Initial values for v_la and v_lv
        z0 = [self.V_la, self.V_lv]
        # Time between t[i-1] and t[i]
        t = [a, b]
        # Solve ode
        z = odeint(self.ODE_blood_volume_Left_Atria_Ventricles, z0, t)
        # Store values after solve ode
        self.V_la = z[1, 0]
        self.V_lv = z[1, 1]

        return (self.V_la, self.V_lv)



    def ODE_pulmonary_circulation(self, x, t):

        # Store the values in x array
        P_pulmonary_artery = x[0]
        P_pulmonary_vein = x[1]
        Q_pulmonary_artery = x[2]
        Q_pulmonary_vein = x[3]


        # Define ode formulation for pulmonary artery pressure
        if self.P_rv >= P_pulmonary_artery:
            der_pulmonary_press_artery = ((self.CQ_puv * sqrt(abs((self.e_v * self.E_rv_A + self.E_rv_B) * (self.V_rv - self.V_rv_0)) - P_pulmonary_artery)) - Q_pulmonary_artery) / self.C_pulmonary_artery
        else:
            der_pulmonary_press_artery = -1 * Q_pulmonary_artery / self.C_pulmonary_artery

        # Define ode formulation for pulmonary vein pressure
        der_pulmonary_press_vein = (Q_pulmonary_artery - Q_pulmonary_vein) / self.C_pulmonary_vein

        # Define ode formulation for pulmonary artery flow
        der_pulmonary_flow_artery = (P_pulmonary_artery - P_pulmonary_vein - Q_pulmonary_artery * self.R_pulmonary_artery) / self.I_pulmonary_artery

        # Define ode formulation for pulmonary vein flow
        der_pulmonary_flow_vein = (P_pulmonary_vein - ((self.e_a * self.E_la_A + self.E_la_B) * (self.V_la - self.V_la_0)) - Q_pulmonary_vein * self.R_pulmonary_vein) / self.I_pulmonary_vein

        return (der_pulmonary_press_artery, der_pulmonary_press_vein, der_pulmonary_flow_artery, der_pulmonary_flow_vein)



    def Solve_ODE_pulmonary_circulation(self, a, b):
        # Initial values for P_pulmonary_artery, P_pulmonary_vein, Q_pulmonary_artery, and Q_pulmonary_vein
        x0 = [self.P_pulmonary_artery, self.P_pulmonary_vein, self.Q_pulmonary_artery, self.Q_pulmonary_vein]
        # Time between t[i-1] and t[i]
        t = [a, b]
        # Solve ode
        x = odeint(self.ODE_pulmonary_circulation, x0, t)
        # Store values after solve ode
        self.P_pulmonary_artery = x[1, 0]
        self.P_pulmonary_vein = x[1, 1]
        self.Q_pulmonary_artery = x[1, 2]
        self.Q_pulmonary_vein = x[1, 3]

        return (self.P_pulmonary_artery, self.P_pulmonary_vein, self.Q_pulmonary_artery, self.Q_pulmonary_vein)




def plot_results(dic_array, x_lable, y_lable, plot_title, dic_array_2={}):
    ''' This function plot the result found solving systems'''

    if dic_array_2 == {}:
        dic_vals = dic_array.items()
        x, y = zip(*dic_vals)
        plt.plot(x, y)
        plt.xlabel(x_lable)
        plt.ylabel(y_lable)
        plt.title(plot_title)
    else:
        dic_vals = dic_array.items()
        x, y = zip(*dic_vals)
        dic_vals_2 = dic_array_2.items()
        x_2, y_2 = zip(*dic_vals_2)
        plt.plot(x, y)
        plt.plot(x_2, y_2)
        plt.xlabel(x_lable)
        plt.ylabel(y_lable)
        plt.title(plot_title)

    # plt.legend()

    plt.show()


def print_results(*arg):
    ''' This function print output values'''
    print("e_a:" + str([(k, arg[0][k]) for k in arg[0]]))
    print("e_v:" + str([(k, arg[1][k]) for k in arg[1]]))

    print("P_ra" + str([(k, arg[2][k]) for k in arg[2]]))
    print("P_rv:" + str([(k, arg[3][k]) for k in arg[3]]))
    print("P_la:" + str([(k, arg[4][k]) for k in arg[4]]))
    print("P_lv:" + str([(k, arg[5][k]) for k in arg[5]]))

    print("Q_la:" + str([(k, arg[6][k]) for k in arg[6]]))
    print("Q_lv:" + str([(k, arg[7][k]) for k in arg[7]]))
    print("Q_ra:" + str([(k, arg[8][k]) for k in arg[8]]))
    print("Q_rv:" + str([(k, arg[9][k]) for k in arg[9]]))

    print("V_ra:" + str([(k, arg[10][k]) for k in arg[10]]))
    print("V_rv:" + str([(k, arg[11][k]) for k in arg[11]]))
    print("V_la:" + str([(k, arg[12][k]) for k in arg[12]]))
    print("V_lv:" + str([(k, arg[13][k]) for k in arg[13]]))

    print("P_pulmonary_artery:" + str([(k, arg[14][k]) for k in arg[14]]))
    print("P_pulmonary_vein:" + str([(k, arg[15][k]) for k in arg[15]]))
    print("Q_pulmonary_artery:" + str([(k, arg[16][k]) for k in arg[16]]))
    print("Q_pulmonary_vein:" + str([(k, arg[17][k]) for k in arg[17]]))



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # define an object regarding heart_lung class
    heart_lung_obj = Heart_Lung()

    # activation_ventricles and activation_atria
    e_v = {}
    e_a = {}

    # blood_press_Atria_ventricles
    P_ra = {}
    P_rv = {}
    P_la = {}
    P_lv = {}

    # blood_flow_atria_ventricles
    Q_ra = {}
    Q_rv = {}
    Q_la = {}
    Q_lv = {}

    # # blood_volume_Atria_ventricles
    V_ra = {}
    V_rv = {}
    V_la = {}
    V_lv = {}

    # pulmonary pressure and flow
    P_pulmonary_artery = {}
    P_pulmonary_vein = {}
    Q_pulmonary_artery = {}
    Q_pulmonary_vein = {}



    # receive an the first and the last element of interval and number of equal parts from user to create an array
    first_int = int(input("please enter the first value of interval :"))
    end_int = int(input("please enter the last value of interval :"))
    num_split = int(input("please enter the number of splitting interval :"))

    t = np.linspace(first_int, end_int, num_split)
    heart_lung_obj.define_t_array(first_int, end_int, num_split)

    # print(flow_file.readline())
    for i in range(0, len(heart_lung_obj.t)):
        if i == 0:  # In this loop we want to run the heart model
            # define time parameter
            mt = heart_lung_obj.t[i] - heart_lung_obj.T * floor(heart_lung_obj.t[i] / heart_lung_obj.T)
            # print(mt)

            # Calculate activation ventricles value regarding its time
            e_v[i] = heart_lung_obj.activation_ventricles(mt)

            # Calculate activation atria value regarding its time
            e_a[i] = heart_lung_obj.activation_atria(mt)

            # Initial value for v_ra, v_rv, v_la, and v_lv
            V_ra[i] = 20.0e-6
            V_rv[i] = 500.0e-6
            V_la[i] = 20.0e-6
            V_lv[i] = 500.0e-6

            # # Calculate blood pressure in atria and ventricles regarding their time
            [P_ra[i], P_rv[i], P_la[i], P_lv[i]] = heart_lung_obj.Blood_Press_Atria_Ventricles(V_ra[i], V_rv[i],
                                                                                               V_la[i], V_lv[i])

            # Initial values for P_pulmonary_artery, P_pulmonary_vein, Q_pulmonary_artery, and Q_pulmonary_vein
            P_pulmonary_artery[i] = 4000
            P_pulmonary_vein[i] = 0
            Q_pulmonary_artery[i] = 0
            Q_pulmonary_vein[i] = 0


            # Calculate blood flow in atria and ventricles regarding their time
            [Q_ra[i], Q_rv[i], Q_la[i], Q_lv[i]] = heart_lung_obj.Blood_Flow_Atria_ventricles()


        else:

            mt = heart_lung_obj.t[i] - heart_lung_obj.T * floor(heart_lung_obj.t[i] / heart_lung_obj.T)
            # print(mt)

            # Calculate activation ventricles value regarding its time
            e_v[i] = heart_lung_obj.activation_ventricles(mt)

            # Calculate activation atria value regarding its time
            e_a[i] = heart_lung_obj.activation_atria(mt)

            # Solve ode blood volume in right atrium and ventricle regarding its time
            [V_ra[i], V_rv[i]] = heart_lung_obj.Solve_ODE_blood_volume_Right_Atria_Ventricles(t[i - 1], t[i])
            heart_lung_obj.V_ra = V_ra[i]
            heart_lung_obj.V_rv = V_rv[i]

            # Solve ode blood volume in left atrium and ventricle regarding its time
            [V_la[i], V_lv[i]] = heart_lung_obj.Solve_ODE_blood_volume_Left_Atria_Ventricles(t[i - 1], t[i])
            heart_lung_obj.V_la = V_la[i]
            heart_lung_obj.V_lv = V_lv[i]


            # Calculate blood pressure in atria and ventricles regarding their time
            [P_ra[i], P_rv[i], P_la[i], P_lv[i]] = heart_lung_obj.Blood_Press_Atria_Ventricles(V_ra[i], V_rv[i],
                                                                                               V_la[i], V_lv[i])

            # Solve ode equations in pulmonary arteries and veins
            [P_pulmonary_artery[i], P_pulmonary_vein[i], Q_pulmonary_vein[i], Q_pulmonary_artery[i]] = heart_lung_obj.Solve_ODE_pulmonary_circulation(t[i - 1], t[i])
            heart_lung_obj.P_pulmonary_artery = P_pulmonary_artery[i]
            heart_lung_obj.P_pulmonary_vein = P_pulmonary_vein[i]
            heart_lung_obj.Q_pulmonary_artery = Q_pulmonary_artery[i]
            heart_lung_obj.Q_pulmonary_vein = Q_pulmonary_vein[i]



            # Calculate blood flow in atria and ventricles regarding their time
            [Q_ra[i], Q_rv[i], Q_la[i], Q_lv[i]] = heart_lung_obj.Blood_Flow_Atria_ventricles()




    plot_results(e_a, 'Time', 'Value', 'Activation Atria')
    plot_results(e_v, 'Time', 'Value', 'Activation Ventricle')

    plot_results(P_ra, 'Time', 'Value', 'Right Atrium Pressure')
    plot_results(P_rv, 'Time', 'Value', 'Right Ventricle Pressure')
    plot_results(P_la, 'Time', 'Value', 'Left Atrium Pressure')
    plot_results(P_lv, 'Time', 'Value', 'Left Ventricle Pressure')


    plot_results(Q_la, 'Time', 'Value', 'Left Atrium Flow')
    plot_results(Q_lv, 'Time', 'Value', 'Left Ventricle Flow')
    plot_results(Q_ra, 'Time', 'Value', 'Right atrium Flow')
    plot_results(Q_rv, 'Time', 'Value', 'Right Ventricle Flow')


    plot_results(V_ra, 'Time', 'Value', 'blood volume in RA')
    plot_results(V_rv, 'Time', 'Value', 'blood volume in RV')
    plot_results(V_la, 'Time', 'Value', 'blood volume in LA')
    plot_results(V_lv, 'Time', 'Value', 'blood volume in LV')


    plot_results(P_pulmonary_artery, 'Time', 'Value', 'Pulmonary Artery Pressure')
    plot_results(P_pulmonary_vein, 'Time', 'Value', 'Pulmonary Vein Pressure')
    plot_results(Q_pulmonary_artery, 'Time', 'Value', 'Pulmonary Artery Flow')
    plot_results(Q_pulmonary_vein, 'Time', 'Value', 'Pulmonary Vein Flow')



    print_results(e_a, e_v, P_ra, P_rv, P_la, P_lv, Q_la, Q_lv, Q_ra, Q_rv,
                  V_ra, V_rv, V_la, V_lv, P_pulmonary_artery, P_pulmonary_vein, Q_pulmonary_artery, Q_pulmonary_vein)
