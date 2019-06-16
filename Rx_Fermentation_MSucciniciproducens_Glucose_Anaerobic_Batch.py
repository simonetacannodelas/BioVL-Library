# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 10:42:13 2019

@author: simoca
"""

from scipy.integrate import odeint 
#Package for plotting 
import math 
#Package for the use of vectors and matrix 
import numpy as np
import array as arr
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import glob
from random import sample
import random
import time


class MSucciniciproducens_anae:
    def __init__(self, Control = False):
        self.mu_max= 1.324
        self.ms= 0.061
        self.i= 1.301
        self.Yaa= 0.999
        self.Yfa= 1.532
        self.Yla= 0.999
        self.Ysa= 1.310
        self.Yx= 0.765
        self.Kd=0.010
        self.Ks= 1.123
        self.Ki= 88.35
        self.P_crit= 17.23
        self.alpha_aa= 0.626
        self.alpha_fa= 0.665
        self.alpha_la= 0
        self.alpha_sa= 1.619
        self.beta_aa= 0.124
        self.beta_fa= 0.105
        self.beta_la= 0.210
        self.beta_sa= 0.355
        #Initial concentrations
        self.X0= 0.1
        self.SA0= 0
        self.AA0= 0
        self.FA0= 0
        self.LA0= 0
        self.G0= 10
        self.V0= 2
        
        self.t_end=10
        self.t_start=0
        
        #parameters for control, default every 1/24 hours:
        self.Control = Control
        self.coolingOn = True
        self.Contamination=False
        self.steps = (self.t_end - self.t_start)*24
        self.T0 = 30
        self.K_p = 2.31e+01
        self.K_i = 3.03e-01
        self.K_d = -3.58e-03
        self.Tset = 30
        self.u_max = 150
        self.u_min = 0

    def rxn(self, C,t , u, fc):
        #when there is no control, k has no effect
        k=1
        #when cooling is off than u = 0
        if self.coolingOn == False:
            u = 0
        if self.Contamination == True:
            fc=np.random.randint(0,10)
            fc=fc/17
        if self.Control == True :
            #Cardinal temperature model with inflection: Salvado et al 2011 "Temperature Adaptation Markedly Determines Evolution within the Genus Saccharomyces"
            #Strain S. cerevisiae PE35 M
            Topt = 30
            Tmax = 45.48
            Tmin = 5.04 
            T = C[5]
            if T < Tmin or T > Tmax:
                 k = 0
            else:
                 D = (T-Tmax)*(T-Tmin)**2
                 E = (Topt-Tmin)*((Topt-Tmin)*(T-Topt)-(Topt-Tmax)*(Topt+Tmin-2*T))
                 k = D/E 
                      
        #number of components
        ##initialize the overall conversion vector
        r=np.zeros((7,1))
        if C[0]< self.P_crit:
            r[0,0]=(self.mu_max*C[5])/(C[5]+self.Ks+((C[5]**2)/self.Ki))*((1-(C[0])/self.P_crit))**self.i
        elif C[0] >= self.P_crit:
            r[0,0]=-self.Kd*C[0]
        else:
            print('la has cagado')
        r[1,0]=self.alpha_sa*r[0,0]+self.beta_sa*C[0]
        r[2,0]=self.alpha_aa*r[0,0]+self.beta_aa*C[0]
        r[3,0]=self.alpha_fa*r[0,0]+self.beta_fa*C[0]
        r[4,0]=self.alpha_la*r[0,0]+self.beta_la*C[0]
        r[5,0]=(1/self.Yx)*r[0,0]+(1/self.Ysa)*r[1,0]+(1/self.Yaa)*r[2,0]+(1/self.Yfa)*r[3,0]+(1/self.Yla)*r[4,0]+self.ms*C[0]

        #Solving the mass balances
        dXdt = r[0,0]
        dSAdt = r[1,0]*fc
        dAAdt = r[2,0]
        dFAdt = r[3,0]
        dLAdt = r[4,0]
        dGdt= -r[5,0]
        dVdt = 0
        if self.Control == True :
            '''
             dHrxn heat produced by cells estimated by yeast heat combustion coeficcient dhc0 = -21.2 kJ/g
             dHrxn = dGdt*V*dhc0(G)-dEdt*V*dhc0(E)-dXdt*V*dhc0(X)
             (when cooling is working)  Q = - dHrxn -W ,
             dT = V[L] * 1000 g/L / 4.1868 [J/gK]*dE [kJ]*1000 J/KJ
             dhc0(EtOH) = -1366.8 kJ/gmol/46 g/gmol [KJ/g]
             dhc0(Glc) = -2805 kJ/gmol/180g/gmol [KJ/g]
             
            ''' 
            #Metabolic heat: [W]=[J/s], dhc0 from book "Bioprocess Engineering Principles" (Pauline M. Doran) : Appendix Table C.8 
            dHrxndt =   dXdt*C[0]*(-21200) #[J/s]  + dGdt*C[4]*(15580)- dEdt*C[4]*(29710) 
            #Shaft work 1 W/L1
            W = 1*C[0] #[J/S] negative because exothermic  
            #Cooling just an initial value (constant cooling to see what happens)
            #dQdt = -0.03*C[4]*(-21200) #[J/S]   
            #velocity of cooling water: u [m3/h] -->controlled by PID       
            
            #Mass flow cooling water
            M=u/3600*1000 #[kg/s]
            #Define Tin = 5 C, Tout=TReactor
            #heat capacity water = 4190 J/kgK
            Tin = 5
            #Estimate water at outlet same as Temp in reactor
            Tout = C[7]
            cpc = 4190
            #Calculate Q from Eq 9.47
            Q=-M*cpc*(Tout-Tin) # J/s    
            #Calculate Temperature change
            dTdt = -1*(dHrxndt - Q + W)/(C[4]*1000*4.1868) #[K/s]
        else: 
            dTdt = 0
        return [dXdt,dSAdt,dAAdt,dFAdt,dLAdt,dGdt,dVdt, dTdt]
                
    def solve(self):
        #solve normal:
        t = np.linspace(self.t_start, self.t_end, self.steps)
        if self.Control == False :
            u = 0
            fc= 1
            C0 = [self.X0, self.SA0, self.AA0, self.FA0, self.LA0, self.G0, self.V0, self.T0]
            C = odeint(self.rxn, C0, t, rtol = 1e-7, mxstep= 500000, args=(u,fc,))
    
        #solve for Control
        else:
            fc=0
            """
            PID Temperature Control:
            """
            # storage for recording values
            C = np.ones([len(t), 7]) 
            C0 = [self.X0, self.SA0, self.AA0, self.FA0, self.LA0, self.G0,self.V0, self.T0]
            C[0] = C0
            self.ctrl_output = np.zeros(len(t)) # controller output
            e = np.zeros(len(t)) # error
            ie = np.zeros(len(t)) # integral of the error
            dpv = np.zeros(len(t)) # derivative of the pv
            P = np.zeros(len(t)) # proportional
            I = np.zeros(len(t)) # integral
            D = np.zeros(len(t)) # derivative
            
            for i in range(len(t)-1):
                #print(t[i])
                #PID control of cooling water
                dt = t[i+1]-t[i]
                #Error
                e[i] = C[i,5] - self.Tset  
                #print(e[i])
                if i >= 1:
                    dpv[i] = (C[i,5]-C[i-1,5])/dt
                    ie[i] = ie[i-1] + e[i]*dt
                P[i]=self.K_p*e[i]
                I[i]=self.K_i*ie[i]
                D[i]=self.K_d*dpv[i]
                
                self.ctrl_output[i]=P[i]+I[i]+D[i]
                u=self.ctrl_output[i]
                if u>self.u_max:
                    u=self.u_max
                    ie[i] = ie[i] - e[i]*dt # anti-reset windup
                if u < self.u_min:
                    u =self.u_min
                    ie[i] = ie[i] - e[i]*dt # anti-reset windup
                #time for solving ODE    
                ts = [t[i],t[i+1]]
                #disturbance
                #if self.t[i] > 5 and self.t[i] < 10:
                #   u = 0                
                #solve ODE from last timepoint to new timepoint with old values              
    
                y =  odeint(self.rxn, C0, ts, rtol = 1e-7, mxstep= 500000, args=(u,fc,))
                #update C0
                C0 = y[-1]
                #merge y to C
                C[i+1]=y[-1]

        return t, C

class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None, width= 9 , height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        fig.subplots_adjust(right = 0.6,top = 0.8)
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        
        FigureCanvas.setSizePolicy(self,
                QSizePolicy.Expanding,
                QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.plot(C ,t)

    def plot(self,t, C):
        time = t
        data = C
        #plot C dependent of t:
        #axe biomass
        ax = self.figure.add_subplot(111)
        #axe succinic acid
        ax2 = ax.twinx()
        #axe acetic acid
        ax3 = ax.twinx()
        #axe formic acid
        ax4 = ax.twinx()
        #axe lactic acid
        ax5=ax.twinx()
        #axe glucose
        ax6=ax.twinx()
        
        ax3.spines['right'].set_position(('axes',1.15))
        ax4.spines['right'].set_position(('axes',1.37))
        ax5.spines['right'].set_position(('axes',1.49))
        ax6.spines['right'].set_position(('axes',1.62))
        
        l1, = ax.plot(time, data[:len(time),0], color = 'g',label = 'Biomass')
        l2, = ax2.plot(time, data[:len(time),1], color = 'r',label = 'Succinic acid')
        l3, = ax3.plot(time, data[:len(time),2]*1000, color = 'y',label = 'Acetic acid')
        l4, = ax4.plot(time, data[:len(time),3]*1000, color = 'y',label = 'Formic acid')
        l5, = ax5.plot(time, data[:len(time),4]*1000, color = 'y',label = 'Lactic acid')        
        l6, = ax6.plot(time, data[:len(time),5], color = 'b',label = 'Glucose')
        
        labels = [l1, l2, l3, l4, l5, l6]
        ax.legend(labels, [l.get_label() for l in labels], loc = 'lower left', bbox_to_anchor = (0,1.05,1,0.5),mode = 'expand', ncol=3)
          
        ax.set_title('')
        ax.set_xlabel('time [h]')
        ax6.set_ylabel('Glucose [g/L]', color = 'b')
        ax2.set_ylabel('Succinic Acid [g/L]', color = 'm')
        ax.set_ylabel('Biomass [g/L]', color = 'g')
        ax3.set_ylabel('Acetic Acid [mg/L]', color = 'c')
        ax4.set_ylabel('Formic Acid [mg/L]', color = 'k')
        ax5.set_ylabel('Lactic Acid [mg/L]', color = 'y')
        
        if hasattr(modelInUse, 't_expfb_start'):
            ax7 = ax.twinx()
            ax7.spines['right'].set_position(('axes',1.57))
            l7, = ax7.plot(time, data[:len(time),6], color = 'grey',label = 'Volume')
            labels = [l1, l2, l3, l4, l5, l6, l7]
            ax.legend(labels, [l.get_label() for l in labels], loc = 'lower left', bbox_to_anchor = (0,1.05,1,0.5),mode = 'expand', ncol=3, fontsize = 'small')
            ax7.set_ylabel('Volume [L]', color = 'grey')
            
            ymin, ymax = ax.get_ylim()
            #seperate phases on the diagramm and plot text which phase is plotted
            ax.text(x = (modelInUse.t_expfb_start/2) , y = ymax + 0.5 ,
                    s= 'Batch phase', fontsize = 7, horizontalalignment='center' )
            ax.axvline(x= modelInUse.t_expfb_start, color = 'black', linestyle = '--')
            ax.text(x = (modelInUse.t_expfb_start+modelInUse.t_constfb_start)/2  , y = ymax + 0.5, 
                        s= 'exp FB', fontsize = 7, horizontalalignment='center')
            
            ax.axvline(x=modelInUse.t_constfb_start, color = 'black', linestyle = '--')
            ax.text(x = (modelInUse.t_constfb_start+modelInUse.t_end)/2 , y = ymax + 0.5,
                    s= 'const FB', fontsize = 7, horizontalalignment='center')
        if hasattr(modelInUse, 't_feed_start'):
            ax7 = ax.twinx()
            ax7.spines['right'].set_position(('axes',1.55))
            l7, = ax7.plot(time, data[:len(time),6], color = 'grey',label = 'Volume')
            labels = [l1, l2, l3, l4, l5, l6, l7]
            ax.legend(labels, [l.get_label() for l in labels], loc = 'lower left', bbox_to_anchor = (0,1.05,1,0.5),mode = 'expand', ncol=3, fontsize = 'small')
            ax5.set_ylabel('Volume [L]', color = 'grey')
            
            
            ymin, ymax = ax.get_ylim()
            ax.text(x = (modelInUse.t_feed_start/2) , y = ymax + 0.5 ,
                    s= 'Batch phase', fontsize = 7, horizontalalignment='center' )
            ax.axvline(x= modelInUse.t_feed_start, color = 'black', linestyle = '--')
            ax.text(x = (modelInUse.t_feed_start+modelInUse.t_end)/2  , y = ymax + 0.5, 
                        s= 'CSTR phase', fontsize = 7, horizontalalignment='center')
                        
        self.draw()

#f=MSucciniciproducens_anae()
#C=f.solve()[1]
#plt.plot(f.solve()[0], C[0],'g-')
#plt.plot(f.solve()[0], C[5],'k-')
#plt.show()
