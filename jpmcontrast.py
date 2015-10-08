# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
from scipy.constants import hbar, h, e, pi

phi0 = h/(2*e)

#Junction parameters

Jc = 35;
A = 500;

Lm = 200;
Cm = 5;

ib = 0.995;

x0  = [Jc, A, ib, Lm, Cm];


def JPMFOM(Ic, Cj, Zc, w, i):
#    Plasma frequency    
    wp = 2**(.25) * np.sqrt(2*pi*Ic/(Cj*phi0)) * (1-i)**(.25)
    
#    Barrier Height
    dU = 4*Ic*phi0 / (3*np.sqrt(2)*pi) * (1-i)**(3./2)
    
    w01 = wp*(1-5*hbar*wp/(36*dU))
    N = dU / (hbar*wp)

    g0 = wp/(2*pi) * np.sqrt(432*dU/(hbar*wp)) * np.exp(-36*dU/(5*hbar*wp))  
    g1 = wp/(2*pi) * (432*dU/(hbar*wp))**(3/2) * np.exp(-36*dU/(5*hbar*wp))/(np.sqrt(pi))
    
    Zj = 1./(wp*Cj)
    gTL = 0.25*np.real(1/Zc)/Cj

    Zi = 1./(1/Zc + 1/50.)
    gin = np.real(1./Zi)/Cj

    Delta = w01 - w
    eta = gTL*(1-g0/(gTL + g1 + gin)) * (g1 + g0)/((.5*(g1+gTL+g0+gin))**2 + Delta**2)    
    return eta, g0
    
#def JPMContrast(x):
    
    
x = x0    
w = 2*pi*5e9; #5GHz
Ic = x[0]*x[1]*1e-8
Cj = x[1]*50e-15
ib = x[2]
t = 10e-9

Zm = 1j*w*x[3]*1e-12 + 1./(1j*w*x[4]*1e-12 * 1./50)

[eta, g0] = JPMFOM(Ic, Cj, Zm, w, ib)

kappa = 16e6
ncrit = 100.
lbda = 2.*kappa*ncrit
C = 1-(1-np.exp(-lbda*eta*t)) * np.exp(-g0*t)

C
    
#    return C
#    
#C = JPMContrast(x0)


#function C = JPMContrast(x)
#    
#    w = 2*pi*5e9;
#    %x = [Jc, Area, I, Lm, Cm]
#    Ic = x(1)*x(2)*1e-8;
#    Cj = x(2)*50e-15;
#    ib = x(3);
#    t  = 10e-9;
#    
#    Zm = 1i*w*x(4)*1e-12 + 1./(1i*w*x(5)*1e-12 + 1/50);
#    
#    [eta, g0] = JPM_FOM(Ic, Cj, Zm, w, ib);
#    kappa = 16e6;
#    ncrit = 100;
#    lambda = 2*kappa*ncrit;
#    C = 1-(1-exp(-lambda*eta*t)).*exp(-g0*t);
#
#end


#function [eta, g0] = JPM_FOM(Ic, Cj, Zc, w, i)
#    %Constants
#    phi0 = 2.067e-15;
#    hbar = 1.0545e-34;
#    %%%%WKB FORMULAE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#    %junction plasma frequency
#    wp = 2^(1/4)*sqrt(2*pi*Ic/Cj/phi0)*(1-i).^(1/4);
#
#    %junction barrier height
#    dU = 4*Ic*phi0/3/sqrt(2)/pi*(1-i).^(3/2);
#    w01 = wp.*(1-5*hbar*wp/36./dU);
#    N = dU./(hbar*wp);
#
#    %tunneling from ground, excited state
#    g0 = wp/2/pi.*sqrt(432*dU./(hbar*wp)).*exp(-36*dU./(5*hbar*wp));
#    g1 = wp/2/pi.*(432*dU./(hbar*wp)).^(3/2).*exp(-36*dU./(5*hbar*wp))/sqrt(pi);
#
#    %%%OTHER JUNCTION PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#    %junction inductance
#    Zj = 1./(wp*Cj);
#    gTL = 0.25*real(1./Zc)/Cj;
#    
#    Zi = 1./(1./Zc + 1/50);
#    gin = real(1./Zi)/Cj;
#    
#    Delta = w01 - w;
#    eta = gTL.*(1 - g0./(gTL+g1+gin)).*(g1 + g0)./((0.5*(g1+gTL+g0+gin)).^2 + Delta.^2);
#    
#    
#end