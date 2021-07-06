import math
import scipy.integrate as integrate
import scipy.special as special
from math import log as ln

########################################
# inputs
print('Enter the di-electric constant:')
er=2.2; #int(input());
print('Enter the substrate thickness (in mm)')
h=1.6; #int(input());
print('Enter the frequency (GHz):')
f=2.4; #int(input());
print('Enter your estimated feed width (in mm):')
w=6; #int(input());


#######################################
f=f*1e9;
u0 = 1.25663706e-6; #permeability of free space
e0 = 8.85418782e-12; #permittivity of free space
eta = 3.767303135e2; #wave impedance of free space
c = 299792458; #speed of light


wid=(c/(math.sqrt((er+1)/2)*2*f))*1000;        #in mm
e_eff=((er+1)/2)+ (((er-1)/2)* (1+((12*h)/wid))**-0.5);
l_eff=(c/(2*f*math.sqrt(e_eff)))*1000;
del_l=(((e_eff+0.3)*((wid/h)+0.264))/((e_eff-0.258)*((wid/h)+0.8)))*(0.412*h);  #in mm
L=l_eff-(2*del_l);
increment = 0.05; # change this increament value for faster iterations


# width and length
lda=c/f*1000;
k=(2*math.pi)/lda
x=k*wid;
sinintx, consintx = special.sici(x);
i1 = -2 + math.cos(x) + (x*sinintx) + (math.sin(x)/x) ;
g1 = i1/(120*math.pi**2);          #Conductance 1
a = lambda th: (((math.sin((x/2)*math.cos(th))/math.cos(th))**2)*(special.jv(0,(k*L*math.sin(th))))*(math.sin(th))**3);
a1 = integrate.quad(a, 0, math.pi)[0]
g12=a1/(120*math.pi**2);     #in siemens #special.jv(0,(k*L*math.sin(th)))Conductance 2
r_in=1/(2*(g1+g12));    #in ohms
inset=(L/math.pi)*(math.acos(math.sqrt(50/r_in)));        #in mm
print("The width is: ", wid, "mm")
print("The length is: ", L, "mm")
print("The inset feed point is: ", inset, "mm")


# impedance width:

def testimpedance(w,e_eff):
    #calculation of effective dielectric constant
    global w_corrected
    er_eff = e_eff;
    k_0 = 2*math.pi*f/c;
    beta = k_0*math.sqrt(er_eff);
    testwh = w/h;
    #calculation of characteristic impedance
    if testwh >= 1:
      z1 = math.log(w/h + 1.444);
      z2 = w/h + 1.393 + 0.667*z1;
      z3 = math.sqrt(er_eff)*z2;
      Z_0 = eta/z3;
    else:
      z1 = math.log(8*h/w + w/(4*h));
      z2 = 60/math.sqrt(er_eff);
      Z_0 = z2*z1;


    if Z_0<49.8:
        #print("iterating ...Reducing the feed width. Current impedance of the feed: {} Ohms" .format(Z_0))
        w = w - increment;
        Z_0 = testimpedance(w,e_eff)[0];
    elif Z_0>50.2:
        #print("iterating...Increamath.sint the feed width. Current impedance of the feed: {} Ohms".format(Z_0))
        w = w + increment;
        Z_0 = testimpedance(w,e_eff)[0];
    else:
        w_corrected = w;
        print("Iterations complete...")

        print("Determined width of the feed (in mm): ", w)
        pass
    return Z_0, w

Z_0, w = testimpedance(w,e_eff)
print("Current impedance of the feed (in mm): ", Z_0)

#inse gap
g = ((c*(4.65e-12))/(math.sqrt(2*e_eff)*(f*(1e-9))));
print("the inset gap", g*1000)
