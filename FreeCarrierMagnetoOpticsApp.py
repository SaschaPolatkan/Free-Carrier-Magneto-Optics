import re #for reading numbers out of a string
import matplotlib
from tkinter import *
from tkinter import filedialog
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np 
#matplotlib.rcParams['text.usetex'] = True
#do not enable tex! Slows down program by a lot. 

np.warnings.filterwarnings('ignore')


wc=0.0 #universal
kl=16.0 #universal
wp=[1.0] #carrier dependent (determines mass)
w0=[0.0] #carrier dependent 
theta=[0.0] #technically dependent on shape of the fermi-surface. Here a spherical one is assumed.

#initial values:
nu_min_exponent=-10
nu_max_exponent=1
nu_res=1000

nu_sliderVar=[-10.0]
nu=[0.0] #carrier dependent (scattering)

class ScaleNew(Scale):
	"""a type of Scale where the left click is hijacked to work like a right click"""
	def __init__(self, master=None, **kwargs):
		Scale.__init__(self, master, **kwargs)
		self.bind('<Button-1>', self.set_value)

	def set_value(self, event):
		self.event_generate('<Button-3>', x=event.x, y=event.y)
		return 'break'

CS=1
NCarriers=1
class DropDownUpdateable():
	def __init__(self, parent):
		self.parent = parent
		self.options = ['1']

		self.om_variable = StringVar(self.parent)
		self.om_variable.set(self.options[0])
		self.om_variable.trace('w', self.option_select)

		self.om = OptionMenu(self.parent, self.om_variable, *self.options)
		self.om.grid(column=62, row=0)
		
		
		self.label = Label(self.parent, text='Number of Carriers')
		self.entry = Entry(self.parent)
		self.entry.bind('<Return>',self.update_option_menu)
		self.button = Button(self.parent, text='Update', command=self.update_option_menu)

		self.label.grid(column=63, row=0)
		self.entry.grid(column=63, row=1)
		self.button.grid(column=63, row=2)


	def update_option_menu(self, *args):
		i=1
		global NCarriers
		global wp
		global w0
		global nu
		global theta
		global nu_sliderVar
		global nu_min_exponent
		try:
			NCarriers=int(self.entry.get())
		except:
			print("")
		
		ListString=[]
		while i<=NCarriers:
			ListString.append(str(i))
			i+=1
		self.options=ListString
		self.entry.delete(0, 'end')
		menu = self.om["menu"]
		menu.delete(0, "end")
		for string in self.options:
			menu.add_command(label=string, command=lambda value=string: self.om_variable.set(value))
			i+=1
        
		wp=[1.0]*NCarriers
		nu=[0.0]*NCarriers
		w0=[0.0]*NCarriers
		theta=[0.0]*NCarriers
		nu_sliderVar=[nu_min_exponent]*NCarriers
		
		

	def option_select(self, *args):
		global CS
		CS=int(self.om_variable.get()) #CarrierSelected
		self.label = Label(self.parent, text='Chosen Carrier: '+str(CS))
		self.label.grid(column=63, row=3)
		update_sliders()
	
	


root = Tk()
DropDownUpdateable(root)


root.wm_title("Parabolic Band")
root.configure(background='white')
fig = plt.Figure()
canvas = FigureCanvasTkAgg(fig, root)
#canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
canvas.get_tk_widget().grid(row=1, column=0, rowspan=40, columnspan = 40)
#function embedds the canvas in the tkinter window. this includes the slider.


wc_label = Text(root, width=9, height=1, borderwidth=0,background=root.cget("background"))
wc_label.tag_configure("subscript", offset=-4)
wc_label.insert("insert", u"\u03c9", "", "c", "subscript",", from:")
wc_label.configure(state="disabled")
wc_label.grid(row=2, column=41, columnspan=2,padx=(0,0), sticky=W)

wc_label2 = Text(root, width=4, height=1, borderwidth=0,background=root.cget("background"))
wc_label2.tag_configure("subscript", offset=-4)
wc_label2.insert("insert", "to:")
wc_label2.configure(state="disabled")
wc_label2.grid(row=2, column=45,columnspan=1,padx=(0,0), sticky=W)

wc_info = Text(root, width=36, height=2, borderwidth=0,background=root.cget("background"))
wc_info.tag_configure("subscript", offset=-4)
wc_info.insert("insert", "in units of "+u"\u03c9", "", "p", "subscript", " of first carrier")
wc_info.configure(state="disabled")
wc_info.grid(row=3, column=41,columnspan=20,padx=(0,0), sticky=W)



#initial values:
wc_min=0
wc_max=5


def update_wclimits():
	global wc_min 
	global wc_max
	wc_smin=wc_entrymin.get()
	if(wc_smin!=""): 
		wc_min=float(wc_smin)	
		print('wc_min:  '+str(wc_min)+' wp')
	wc_smax=wc_entrymax.get()
	if(wc_smax!=""): 
		wc_max=float(wc_smax)	
		print('wc_max:  '+str(wc_max)+' wp')
	wc_slider.configure(from_=wc_min, to=wc_max, resolution=abs(wc_max-wc_min)/200)
		
wc_entrymin = Entry(root, width=5, relief=FLAT, highlightbackground="#fb0", highlightthickness=2)
wc_entrymin.grid(row =2, column = 43,columnspan=2,padx=(0,40), sticky=W)
wc_entrymin.bind("<Return>", (lambda event: update_wclimits()))

wc_entrymax = Entry(root, width=5, relief=FLAT, highlightbackground="#fb0", highlightthickness=2)
wc_entrymax.grid(row =2, column = 46,padx=(0,0), sticky=W)
wc_entrymax.bind("<Return>", (lambda event: update_wclimits()))




theta_label = Text(root, width=9, height=1, borderwidth=0,background=root.cget("background"))
theta_label.tag_configure("subscript", offset=-4)
theta_label.insert("insert", u"\u03b8, from:")
theta_label.configure(state="disabled")
theta_label.grid(row=6, column=41,columnspan=2,padx=(0,0), sticky=W)

theta_label = Text(root, width=4, height=1, borderwidth=0,background=root.cget("background"))
theta_label.tag_configure("subscript", offset=-4)
theta_label.insert("insert", "to:")
theta_label.configure(state="disabled")
theta_label.grid(row=6, column=45,columnspan=1,padx=(0,0), sticky=W)

theta_info = Text(root, width=28, height=2, borderwidth=0,background=root.cget("background"))
theta_info.tag_configure("subscript", offset=-4)
theta_info.insert("insert", "in degrees. First carrier\ndefines angle relative to B.")
theta_info.configure(state="disabled")
theta_info.grid(row=7, column=41,columnspan=6,padx=(0,0), sticky=W)

theta_min=0
theta_max=90

def update_thetalimits():
	global theta_min 
	global theta_max
	theta_smin=theta_entrymin.get()
	if(theta_smin!=""): 
		theta_min=float(theta_smin)	
		print('theta_min:  '+str(theta_min))
	theta_smax=theta_entrymax.get()
	if(theta_smax!=""): 
		theta_max=float(theta_smax)	
		print('theta_max:  '+str(theta_max))
	theta_slider.configure(from_=theta_min, to=theta_max, resolution=abs(theta_max-theta_min)/180)

theta_entrymin = Entry(root, width=5, relief=FLAT, highlightbackground="#fb0", highlightthickness=2)
theta_entrymin.grid(row =6, column = 43,columnspan=2,padx=(0,40), sticky=W)
theta_entrymin.bind("<Return>", (lambda event: update_thetalimits()))

theta_entrymax = Entry(root, width=5, relief=FLAT, highlightbackground="#fb0", highlightthickness=2)
theta_entrymax.grid(row =6, column = 46,padx=(0,0), sticky=W)
theta_entrymax.bind("<Return>", (lambda event: update_thetalimits()))


thetalock=BooleanVar() 
thetalock.set(True)
thetalock_checkbox=Checkbutton(root, text="Lock Theta", variable=thetalock, bg='white').grid(row =8, column = 41, sticky=W)




kl_label = Text(root, width=9, height=1, borderwidth=0,background=root.cget("background"))
kl_label.tag_configure("subscript", offset=-4)
kl_label.insert("insert", "kl, from:")
kl_label.configure(state="disabled")
kl_label.grid(row=10, column=41,columnspan=2,padx=(0,0), sticky=W)

kl_label = Text(root, width=4, height=1, borderwidth=0,background=root.cget("background"))
kl_label.tag_configure("subscript", offset=-4)
kl_label.insert("insert", "to:")
kl_label.configure(state="disabled")
kl_label.grid(row=10, column=45,columnspan=1,padx=(0,0), sticky=W)



kl_min=0
kl_max=200

def update_kllimits():
	global kl_min 
	global kl_max
	kl_smin=kl_entrymin.get()
	if(kl_smin!=""): 
		kl_min=float(kl_smin)	
		print('kl_min:  '+str(kl_min))
	kl_smax=kl_entrymax.get()
	if(kl_smax!=""): 
		kl_max=float(kl_smax)	
		print('kl_max:  '+str(kl_max))
	kl_slider.configure(from_=kl_min, to=kl_max, resolution=abs(kl_max-kl_min)/200)

kl_entrymin = Entry(root, width=5, relief=FLAT, highlightbackground="#fb0", highlightthickness=2)
kl_entrymin.grid(row =10, column = 43,columnspan=2,padx=(0,40), sticky=W)
kl_entrymin.bind("<Return>", (lambda event: update_kllimits()))

kl_entrymax = Entry(root, width=5, relief=FLAT, highlightbackground="#fb0", highlightthickness=2)
kl_entrymax.grid(row =10, column = 46,padx=(0,0), sticky=W)
kl_entrymax.bind("<Return>", (lambda event: update_kllimits()))
kl=16




nu_label = Text(root, width=9, height=1, borderwidth=0,background=root.cget("background"))
nu_label.tag_configure("subscript", offset=-4)
nu_label.insert("insert", u"\u03bd, from:")
nu_label.configure(state="disabled")
nu_label.grid(row=6, column=62,columnspan=2,padx=(0,0), sticky=W)

nu_label = Text(root, width=4, height=1, borderwidth=0,background=root.cget("background"))
nu_label.tag_configure("subscript", offset=-4)
nu_label.insert("insert", "to:")
nu_label.configure(state="disabled")
nu_label.grid(row=6, column=65,columnspan=1,padx=(0,0), sticky=W)

nu_info = Text(root, width=14, height=2, borderwidth=0,background=root.cget("background"))
nu_info.tag_configure("subscript", offset=-4)
nu_info.configure(state="disabled")
nu_info.grid(row=7, column=62,columnspan=3,padx=(0,0), sticky=W)


def update_nulimits():
	global nu_res 
	global nu_min_exponent 
	global nu_max_exponent
	nu_smin=nu_entrymin.get()
	if(nu_smin!=""): 
		nu_min_exponent=float(nu_smin)	
		print('nu_min_exponent:  '+str(nu_min_exponent))
	nu_smax=nu_entrymax.get()
	if(nu_smax!=""): 
		nu_max_exponent=float(nu_smax)	
		print('nu_max_exponent:  '+str(nu_max_exponent))
	nu_slider.configure(from_=nu_min_exponent, to=nu_max_exponent, resolution=abs(nu_max_exponent-nu_min_exponent)/nu_res)

nu_entrymin = Entry(root, width=5, relief=FLAT, highlightbackground="#fb0", highlightthickness=2)
nu_entrymin.grid(row =6, column = 63,columnspan=1,padx=(40,30), sticky=W)
nu_entrymin.bind("<Return>", (lambda event: update_nulimits()))

nu_entrymax = Entry(root, width=5, relief=FLAT, highlightbackground="#fb0", highlightthickness=2)
nu_entrymax.grid(row =6, column = 66,padx=(0,0), sticky=W)
nu_entrymax.bind("<Return>", (lambda event: update_nulimits()))



wp_label = Text(root, width=9, height=1, borderwidth=0,background=root.cget("background"))
wp_label.tag_configure("subscript", offset=-4)
wp_label.insert("insert", u"\u03c9", "", "p", "subscript",", from:")
wp_label.configure(state="disabled")
wp_label.grid(row=10, column=62, columnspan=2,padx=(0,0), pady=(0,0), sticky=W)

wp_label2 = Text(root, width=4, height=1, borderwidth=0,background=root.cget("background"))
wp_label2.tag_configure("subscript", offset=-4)
wp_label2.insert("insert", "to:")
wp_label2.configure(state="disabled")
wp_label2.grid(row=10, column=65,columnspan=1,padx=(0,0), sticky=W)

wp_info = Text(root, width=36, height=2, borderwidth=0,background=root.cget("background"))
wp_info.tag_configure("subscript", offset=-4)
wp_info.insert("insert", "in units of "+u"\u03c9", "", "p", "subscript", " of first carrier")
wp_info.configure(state="disabled")
wp_info.grid(row=11, column=62,columnspan=20,padx=(0,0), sticky=W)



#initial values:
wp_min=0
wp_max=5



def update_wplimits():
	global wp_min 
	global wp_max
	wp_smin=wp_entrymin.get()
	if(wp_smin!=""): 
		wp_min=float(wp_smin)	
		print('wp_min:  '+str(wp_min)+' wp')
	wp_smax=wp_entrymax.get()
	if(wp_smax!=""): 
		wp_max=float(wp_smax)	
		print('wp_max:  '+str(wp_max)+' wp1')
	wp_slider.configure(from_=wp_min, to=wp_max, resolution=abs(wp_max-wp_min)/200)
		
wp_entrymin = Entry(root, width=5, relief=FLAT, highlightbackground="#fb0", highlightthickness=2)
wp_entrymin.grid(row =10, column = 63,columnspan=1,padx=(40,30), sticky=W)
wp_entrymin.bind("<Return>", (lambda event: update_wplimits()))

wp_entrymax = Entry(root, width=5, relief=FLAT, highlightbackground="#fb0", highlightthickness=2)
wp_entrymax.grid(row =10, column = 66,padx=(0,0), sticky=W)
wp_entrymax.bind("<Return>", (lambda event: update_wplimits()))





def update_wlimits():
	global w_min 
	global w_max
	w_smin=w_entrymin.get()
	if(w_smin!=""): 
		w_min=float(w_smin)	
		print('w_min:  '+str(w_min))
	w_smax=w_entrymax.get()
	if(w_smax!=""): 
		w_max=float(w_smax)	
		print('w_max:  '+str(w_max))
	ax.axis([w_min,w_max,params_float[2]*10**z,params_float[3]*10**z])
	fig.canvas.draw_idle()
	


	




ax=fig.add_subplot(111)
title = Label(root,  background='white', font=20)
title.grid(row =0, column=0, columnspan=40)


fig.subplots_adjust(bottom=0.25, top = 0.98)




#write function that shows left and right handedness as colorcode. near to wp they swap...encode ellipticity as color! around wp they're not perfectly circular anymore!!
#choose between: plot both circular pols, or just one of each, where only one sign in Im(Ex/Ey) is plotted (but still ellipticity accounted for!)

#initial params_float:

	

parameters = Entry(root, width=21, relief=FLAT, highlightbackground="blue", highlightthickness=2)
parameters.insert(0, "0, 6, 0, 8")
parameters.grid(row=38, column=41,columnspan=5,padx=(0,0), sticky=W)
param_label1= Label(root, text='Adjust XY-Range', background='white', font=("Courier", 10))
param_label1.grid(row=37, column=41,columnspan=5,padx=(0,0), sticky=W)
param_info= Label(root, text='Xmin,Xmax,Ymin,Ymax', background='white', font=("Courier", 8))
param_info.grid(row=39, column=41,columnspan=5,padx=(0,0), sticky=W)


phi=np.linspace(0,2.*np.pi, 100)



w_min=0
w_max=5
Rmin=0
Rmax=1 
N=1000
params_float=[w_min, w_max, -1, 2] #initial values
params_float_k=[0, 20, -40, 40] #not used yet. the k1, k2 values have a larger range than R. thus this range could be used when plotting them.

params_float_zoomed=params_float

getcoord=[0.,0.]


CPButtonCheck=0
FButtonCheck=0
R1ButtonCheck=0
R2ButtonCheck=0
RtotButtonCheck=0
RtotrelButtonCheck=0
k1ButtonCheck=0
k2ButtonCheck=0
ExOverEy1ButtonCheck=0
ExOverEy2ButtonCheck=0
CPsave=0
Fsave=0 #FaradayPerLength
R1save=0
R2save=0
Rtotsave=0
Rtotrelsave=0
k1save=0
k2save=0
ExOverEy1save=0
ExOverEy2save=0

plotFcondition=0
plotR1condition=0
plotR2condition=0
plotRtotcondition=0
plotRtotrelcondition=0
plotk1condition=0
plotk2condition=0
plotExOverEy1condition=0
plotExOverEy2condition=0
plotCPcondition=0

cc=0 #carrier counter, counting variable for going through the lists of carrier-dependent variables

def replot(w_min, w_max, N):
	global plotFcondition
	global plotR1condition
	global plotR2condition
	global plotRtotcondition
	global plotRtotrelcondition
	global plotk1condition
	global plotk2condition
	global plotExOverEy1condition
	global plotExOverEy2condition
	global plotCPcondition
	global lineF
	global lineR1
	global lineR2
	global lineRtot
	global lineRtotrel
	global linek1
	global linek2
	global lineExOverEy1
	global lineExOverEy2
	global lineCP1
	global lineCP2
	
	global getcoord
	global phi
	
	
	global CPButtonCheck
	global FButtonCheck
	global R1ButtonCheck
	global R2ButtonCheck
	global RtotButtonCheck
	global RtotrelButtonCheck
	global k1ButtonCheck
	global k2ButtonCheck
	global ExOverEy1ButtonCheck
	global ExOverEy2ButtonCheck
	global CPsave
	global Fsave #FaradayPerLength
	global R1save
	global R2save
	global Rtotsave
	global Rtotrelsave
	global k1save
	global k2save
	global ExOverEy1save
	global ExOverEy2save
	global wp 
	global w0 
	global nu
	global theta
	global thetalock
	global cc 
	global CS
	
	global w
	global R1
	global R2 
	global Rtot
	global ExOverEy1
	global ExOverEy1imag
	global ExOverEy2
	global ExOverEy2imag
	global FaradayPerLength
	
	if len(wp)!=NCarriers:
		wp=[]
		nu=[]
		i_wp=0
		while i_wp<=NCarriers:
			wp.append(1.0)
			nu.append(0.0)
	
	w= np.linspace(w_min, w_max, N)
	
	
	#swapper=1
	#when swapper is set to 1, then R1 = R1 and R2 = R2 
	#when swapper is set to (wp-w)/abs(wp-w), then R1 = R+ and R2 = R-
	#note, at an angle theta>0 there is still an unavoidable singularity at wp. 
	#this throws the definitions of R+ and R- off around w=wp. check Im(Ex/Ey) to understand what is going on.
	k1=kl
	k2=kl
	k0=kl
	kxx=kl
	kyy=kl
	kzz=kl
	kxy=0
	kxz=0
	kyz=0
	
	
	cc=0
	
	
	while cc<=NCarriers-1:
		mrel=wp[0]**2/wp[cc]**2 #=mc/m1**2 since wp ~ 1/sqrt(m)
		swapper=(wp[cc]-w)/abs(wp[cc]-w) #this allows to swap between k1 and k2 definitions at wp, thus giving R+ and R- rather than R1 and R2, in regards to Paliks definitions 
		#see page 1230 in palik/furdyna:
		k1+=-kl*(2*wp[cc]**2*(w*(w+1.j*nu[cc])-wp[cc]**2))/(w*(2*(w+1.j*nu[cc])*(w*(w+1.j*nu[cc])-wp[cc]**2)-w*(wc*mrel)**2*np.sin(theta[cc]/360*2*np.pi)**2+swapper*(wc*mrel)*(w**2*(wc*mrel)**2*np.sin(theta[cc]/360*2*np.pi)**4+4*(w*(w+1.j*nu[cc])-wp[cc]**2)**2*np.cos(theta[cc]/360*2*np.pi)**2)**0.5))
		k2+=-kl*(2*wp[cc]**2*(w*(w+1.j*nu[cc])-wp[cc]**2))/(w*(2*(w+1.j*nu[cc])*(w*(w+1.j*nu[cc])-wp[cc]**2)-w*(wc*mrel)**2*np.sin(theta[cc]/360*2*np.pi)**2-swapper*(wc*mrel)*(w**2*(wc*mrel)**2*np.sin(theta[cc]/360*2*np.pi)**4+4*(w*(w+1.j*nu[cc])-wp[cc]**2)**2*np.cos(theta[cc]/360*2*np.pi)**2)**0.5))
		k0+=-kl*(2*wp[cc]**2*(w*(w+1.j*nu[cc])-wp[cc]**2))/(w*(2*(w+1.j*nu[cc])*(w*(w+1.j*nu[cc])-wp[cc]**2)))
		#see page 1205/1206 in palik/furdyna::
		
		#primed coordinates, eq (2.43):
		kxx+=kl*(1.0j/w)*wp[cc]**2*((nu[cc]-1.0j*w)**2+(wc*mrel)**2*np.sin(theta[cc])**2)/((nu[cc]-1.0j*w)*((nu[cc]-1.0j*w)**2+(wc*mrel)**2))
		kzz+=kl*(1.0j/w)*wp[cc]**2*((nu[cc]-1.0j*w)**2+(wc*mrel)**2*np.cos(theta[cc])**2)/((nu[cc]-1.0j*w)*((nu[cc]-1.0j*w)**2+(wc*mrel)**2))
		kxy+=kl*(1.0j/w)*wp[cc]**2*((wc*mrel)*np.cos(theta[cc]))/((nu[cc]-1.0j*w)**2+(wc*mrel)**2)
		
		
		cc+=1
	

	k1real=k1.real
	k2real=k2.real
	R1=(abs((1-k1**0.5)/(1+k1**0.5))**2).real
	R2=(abs((1-k2**0.5)/(1+k2**0.5))**2).real
	Rtot=0.5*(R1+R2)
	
	#carriers at arbitrary angles in terms of primed kappa. kappa contains angle dependencies, which makes it suitable for multiple carriers with different fermi surfaces
	ExOverEy1=(kxx-k1)/kxy #k1=N1**2
	ExOverEy2=(kxx-k2)/kxy #k2=N2**2
	
	
	
	#corresponds to equation (3.8a). Note that for measurements along q, with B tipped away from q, the text says that formula (3.18) is to be used. But (3.18) is expressed with unprimed kappas, thus the angle dependency occurs explicitly. 
	#if the primed kappa are used, they contain the angle dependency, and for the ratio Ex'/Ey' we use again equation (3.8a), but with primed kappas instead.

	ExOverEy1imag=ExOverEy1.imag
	ExOverEy2imag=ExOverEy2.imag

	c0=1
	alpha1=w/c0*((abs(k1)-k1real)/2)**0.5
	alpha2=w/c0*((abs(k2)-k2real)/2)**0.5
	FaradayPerLength=0.5*(alpha2-alpha1) 
	#see Palik, Furdyna, eq (4.2) and (4.30)
	#FaradayPerLength=FaradayPerLength*abs(wp-abs(w))/(abs(w)-wp) #due do definitions, there is a discontinuity at wp. this smoothes it out. on both sides, thus abs(w)!
	
	ax.axis(params_float_zoomed)
	
	x_ellipse1=1*np.cos(phi)
	x_ellipse2=1*np.cos(phi)
	y_ellipse1=1*np.sin(phi)
	y_ellipse2=1*np.sin(phi)
	
	fig.canvas.draw_idle()
	#Note that I put c=1!!
	
	if(CPButtonCheck==1):
		canvas.get_tk_widget().bind("<Motion>", Plot_CP_ellipse)
				
	
	if(FButtonCheck==1):
		#print('wc:' + str(wc))
		#FaradayPerLength=0.5*(alpha2-alpha1) #see Palik, Furdyna, eq (4.2) and (4.30)
		ax.axis(params_float_zoomed)
		if(plotFcondition==0):
			lineF, = ax.plot(w, FaradayPerLength) #comma belongs there for some reason, because ax.plot() 'returns a tuple of line objects' ... ??
			plotFcondition=1
		if(plotFcondition>0):
			lineF.set_ydata(FaradayPerLength)
			lineF.set_xdata(w)
		fig.canvas.draw_idle()
	
	if(R1ButtonCheck==1):
		#R1=(abs((1-k1**0.5)/(1+k1**0.5))**2).real #Note that Rplus refers to the choice of the sign in kappa
		#palik formula (4.15) -> R defined as absolute squared of r 
		ax.axis(params_float_zoomed)
		if(plotR1condition==0):
			lineR1, = ax.plot(w, R1) 
			plotR1condition=1
		if(plotR1condition>0):
			lineR1.set_ydata(R1)
			lineR1.set_xdata(w)
		fig.canvas.draw_idle()
		
	if(R2ButtonCheck==1):
		#R1=(abs((1-k1**0.5)/(1+k1**0.5))**2).real #Note that Rplus refers to the choice of the sign in kappa
		#R2=(abs((1-k2**0.5)/(1+k2**0.5))**2).real
		#palik formula (4.15) -> R defined as absolute squared of r 
		ax.axis(params_float_zoomed)
		if(plotR2condition==0):
			lineR2, = ax.plot(w, R2) 
			plotR2condition=1
		if(plotR2condition>0):
			lineR2.set_ydata(R2)
			lineR2.set_xdata(w)
		fig.canvas.draw_idle()
	
	if(RtotButtonCheck==1):
		
		R1=(abs((1-k1**0.5)/(1+k1**0.5))**2).real #Note that Rplus refers to the choice of the sign in kappa
		R2=(abs((1-k2**0.5)/(1+k2**0.5))**2).real
		Rtot=(R1+R2)/2.0
		#palik formula (4.15) -> R defined as absolute squared of r 
		ax.axis(params_float_zoomed)
		if(plotRtotcondition==0):
			lineRtot, = ax.plot(w, Rtot) 
			plotRtotcondition=1
		if(plotRtotcondition>0):
			lineRtot.set_ydata(Rtot)
			lineRtot.set_xdata(w)
		fig.canvas.draw_idle()
	
	if(RtotrelButtonCheck==1):
		R1=(abs((1-k1**0.5)/(1+k1**0.5))**2).real 
		R2=(abs((1-k2**0.5)/(1+k2**0.5))**2).real
		R0=(abs((1-k0**0.5)/(1+k0**0.5))**2).real 
		
		Rtotrel=(R1+R2)/2.0/R0
		
		ax.axis(params_float_zoomed)
		if(plotRtotrelcondition==0):
			lineRtotrel, = ax.plot(w, Rtotrel) 
			plotRtotrelcondition=1
		if(plotRtotrelcondition>0):
			lineRtotrel.set_ydata(Rtotrel)
			lineRtotrel.set_xdata(w)
		fig.canvas.draw_idle()
			
	if(k1ButtonCheck==1):
		ax.axis(params_float_zoomed)
		if(plotk1condition==0):
			linek1, = ax.plot(w, k1real) 
			plotk1condition=1
		if(plotk1condition>0):
			linek1.set_ydata(k1real)
			linek1.set_xdata(w)
		fig.canvas.draw_idle()
		
	if(k2ButtonCheck==1):
		ax.axis(params_float_zoomed)
		if(plotk2condition==0):
			linek2, = ax.plot(w, k2real) 
			plotk2condition=1
		if(plotk2condition>0):
			linek2.set_ydata(k2real)
			linek2.set_xdata(w)
		fig.canvas.draw_idle()
			
	if(ExOverEy1ButtonCheck==1):
		ax.axis(params_float_zoomed)
		if(plotExOverEy1condition==0):
			lineExOverEy1, = ax.plot(w, ExOverEy1imag) 
			plotExOverEy1condition=1
		if(plotExOverEy1condition>0):
			lineExOverEy1.set_ydata(ExOverEy1imag)
			lineExOverEy1.set_xdata(w)
		fig.canvas.draw_idle()
		
	if(ExOverEy2ButtonCheck==1):
		ax.axis(params_float_zoomed)
		if(plotExOverEy2condition==0):
			lineExOverEy2, = ax.plot(w, ExOverEy2imag) 
			plotExOverEy2condition=1
		if(plotExOverEy2condition>0):
			lineExOverEy2.set_ydata(ExOverEy2imag)
			lineExOverEy2.set_xdata(w)
		fig.canvas.draw_idle()
			
	if(FButtonCheck-Fsave<0): #FaradayPerLength
		FButtonCheck=0
		Fsave=0
		lineF.set_ydata([])
		lineF.set_xdata([])
		fig.canvas.draw_idle()
		
	if(R1ButtonCheck-R1save<0):
		R1ButtonCheck=0
		R1save=0
		lineR1.set_ydata([])
		lineR1.set_xdata([])
		fig.canvas.draw_idle()

	if(R2ButtonCheck-R2save<0):
		R2ButtonCheck=0
		R2save=0
		lineR2.set_ydata([])
		lineR2.set_xdata([])
		fig.canvas.draw_idle()
		
	if(RtotButtonCheck-Rtotsave<0):
		RtotButtonCheck=0
		Rtotsave=0
		lineRtot.set_ydata([])
		lineRtot.set_xdata([])
		fig.canvas.draw_idle()
		
	if(RtotrelButtonCheck-Rtotrelsave<0):
		RtotrelButtonCheck=0
		Rtotrelsave=0
		lineRtotrel.set_ydata([])
		lineRtotrel.set_xdata([])
		fig.canvas.draw_idle()
	
	if(k1ButtonCheck-k1save<0):
		k1ButtonCheck=0
		k1save=0
		linek1.set_ydata([])
		linek1.set_xdata([])
		fig.canvas.draw_idle()

	if(k2ButtonCheck-k2save<0):
		k2ButtonCheck=0
		k2save=0
		linek2.set_ydata([])
		linek2.set_xdata([])
		fig.canvas.draw_idle()
		
	if(ExOverEy1ButtonCheck-ExOverEy1save<0):
		ExOverEy1ButtonCheck=0
		ExOverEy1save=0
		lineExOverEy1.set_ydata([])
		lineExOverEy1.set_xdata([])
		fig.canvas.draw_idle()

	if(ExOverEy2ButtonCheck-ExOverEy2save<0):
		ExOverEy2ButtonCheck=0
		ExOverEy2save=0
		lineExOverEy2.set_ydata([])
		lineExOverEy2.set_xdata([])
		fig.canvas.draw_idle()
		
	if(CPButtonCheck-CPsave<0):
		CPButtonCheck=0
		CPsave=0
		lineCP1.set_ydata([])
		lineCP2.set_ydata([])
		lineCP1.set_xdata([])
		lineCP2.set_xdata([])
		fig.canvas.draw_idle()
	
	

	

#Button to plot CircPol, meaning the ratio E_x/E_y in units of i. When 1 -> CircPol, when 0 or infinity (depending on whether E_x/E_y or E_y/E_x is calculated) then linear
def CPButton():
	global CPButtonCheck
	global CPsave
	CPsave=CPButtonCheck
	CPButtonCheck=abs(CPButtonCheck -1) #alternates between 1 and 0 everytime function is called 
	replot(params_float_zoomed[0],params_float_zoomed[1],N)

CPButton=Button(root, text='Pol.', width=5,height=1,relief=FLAT, background="#ccffcc",command=CPButton)
CPButton.grid(row=31, column=46, rowspan=1,padx=(0,0), sticky=W+E)


def FButton():
	global FButtonCheck
	global Fsave #FaradayPerLength
	Fsave=FButtonCheck
	FButtonCheck=abs(FButtonCheck -1) #alternates between 1 and 0 everytime function is called 
	replot(params_float_zoomed[0],params_float_zoomed[1],N)
	
FButton=Button(root, text=u'\u03b8', width=5,height=1,relief=FLAT, background="#ccffcc",command=FButton)
FButton.grid(row=32, column=46, rowspan=1,padx=(0,0), sticky=W+E)

def R1Button():
	global R1ButtonCheck
	global R1save
	R1save=R1ButtonCheck
	R1ButtonCheck=abs(R1ButtonCheck -1) #alternates between 1 and 0 everytime function is called 
	replot(params_float_zoomed[0],params_float_zoomed[1],N)
	
def R2Button():
	global R2ButtonCheck
	global R2save
	R2save = R2ButtonCheck
	R2ButtonCheck=abs(R2ButtonCheck -1) #alternates between 1 and 0 everytime function is called 
	replot(params_float_zoomed[0],params_float_zoomed[1],N)
	
def RtotButton():
	global RtotButtonCheck
	global Rtotsave
	Rtotsave = RtotButtonCheck
	RtotButtonCheck=abs(RtotButtonCheck -1) #alternates between 1 and 0 everytime function is called 
	replot(params_float_zoomed[0],params_float_zoomed[1],N)
	
def RtotrelButton():
	global RtotrelButtonCheck
	global Rtotrelsave
	Rtotrelsave = RtotrelButtonCheck
	RtotrelButtonCheck=abs(RtotrelButtonCheck -1) #alternates between 1 and 0 everytime function is called 
	replot(params_float_zoomed[0],params_float_zoomed[1],N)

RButton=Button(root, text='R1', width=5,relief=FLAT, background="#ccffcc",command=R1Button)
RButton.grid(row=31, column=42, columnspan=1,padx=(0,0), sticky=W+E)
RButton=Button(root, text='R2', width=5,relief=FLAT, background="#ccffcc",command=R2Button)
RButton.grid(row=32, column=42, columnspan=1,padx=(0,0), sticky=W+E)
RButton=Button(root, text='Rtot', width=5,relief=FLAT, background="#ccffcc",command=RtotButton)
RButton.grid(row=31, column=41, columnspan=1, rowspan=1,padx=(0,0), sticky=W+E)
RButton=Button(root, text='Rtot rel.', width=5,relief=FLAT, background="#ccffcc",command=RtotrelButton)
RButton.grid(row=32, column=41, columnspan=1, rowspan=1,padx=(0,0), sticky=W+E)


def ExOverEy1Button():
	global ExOverEy1ButtonCheck
	global ExOverEy1save
	ExOverEy1save=ExOverEy1ButtonCheck
	ExOverEy1ButtonCheck=abs(ExOverEy1ButtonCheck -1)
	replot(params_float_zoomed[0],params_float_zoomed[1],N)

def ExOverEy2Button():
	global ExOverEy2ButtonCheck
	global ExOverEy2save
	ExOverEy2save=ExOverEy2ButtonCheck
	ExOverEy2ButtonCheck=abs(ExOverEy2ButtonCheck -1)
	replot(params_float_zoomed[0],params_float_zoomed[1],N)

ExOverEy1Button=Button(root, text='Im(Ex/Ey)1', width=7,relief=FLAT, background="#ffccb3",command=ExOverEy1Button)
ExOverEy1Button.grid(row=31, column=43,  columnspan=1,padx=(0,0), sticky=W+E)
ExOverEy2Button=Button(root, text='Im(Ex/Ey)2', width=7,relief=FLAT, background="#ffccb3",command=ExOverEy2Button)
ExOverEy2Button.grid(row=32, column=43,  columnspan=1,padx=(0,0), sticky=W+E)


def k1Button():
	global k1ButtonCheck
	global k1save
	k1save=k1ButtonCheck
	k1ButtonCheck=abs(k1ButtonCheck -1)
	replot(params_float_zoomed[0],params_float_zoomed[1],N)

def k2Button():
	global k2ButtonCheck
	global k2save
	k2save=k2ButtonCheck
	k2ButtonCheck=abs(k2ButtonCheck -1)
	replot(params_float_zoomed[0],params_float_zoomed[1],N)

kButton=Button(root, text=u'\u03ba 1', width=6,relief=FLAT, background="#ccffcc",command=k1Button)
kButton.grid(row=31, column=44, columnspan=2,padx=(0,0), sticky=W)
kButton=Button(root, text=u'\u03ba 2', width=6,relief=FLAT, background="#ccffcc",command=k2Button)
kButton.grid(row=32, column=44, columnspan=2,padx=(0,0), sticky=W)



#SLIDERS BELOW PLOT, BECAUSE FUNCTIONS WILL UPDATE PLOT


def exportButton():
	global wc
	global theta
	global nu
	
	global w
	global Rtot
	global R1
	global R2
	global ExOverEy1imag
	global ExOverEy2imag
	file=filedialog.asksaveasfile(mode='w', defaultextension=".txt")
	iterator=0
	if file:
		file.write("wc="+str(wc)+"\t"+"theta="+str(theta)+"\t"+"nu="+str(nu)+"\n")
		file.write("Freq"+"\t"+"Rtot"+"\t"+"R1"+"\t"+"R2"+"\t"+"Im(Ex1/Ey1)"+"\t"+"Im(Ex2/Ey2)"+"\t"+"FaradayPerLength"+"\n")
		while iterator <len(w):
			file.write("%.6f"%w[iterator]+"\t"+"%.6f"%Rtot[iterator]+"\t"+"%.6f"%R1[iterator]+"\t"+"%.6f"%R2[iterator]+"\t"+"%.6f"%ExOverEy1imag[iterator]+"\t"+"%.6f"%ExOverEy2imag[iterator]+"\t"+"%.6f"%FaradayPerLength[iterator]+"\n")
			iterator+=1
		#print(iterator)
		file.close()
	#print(iterator)
	
exportButton=Button(root, text='Export Data', width=11,relief=FLAT, background="#a6c6ff",command=exportButton)
exportButton.grid(row=38, column=46, columnspan=3,padx=(0,0), sticky=W)


def update_wc(var):
	global wc
	global wp
	wc=wp[0]*float(var) 
	replot(params_float_zoomed[0],params_float_zoomed[1],N)
wc_slider= ScaleNew(root, from_=wc_min, to=wc_max, orient=HORIZONTAL, background='white',highlightbackground='white',  sliderlength=5, width=10, length=258, resolution=abs(wc_max-wc_min)/200, sliderrelief=FLAT, troughcolor="#ffddcc", command=update_wc)
wc_slider.grid(row =1, column=41, columnspan=20, padx=(0,30), sticky=W)	

def update_theta(var):
	global theta
	global thetalock
	if(thetalock.get()==True):
		theta_locked=float(var)
		theta=[theta_locked]*NCarriers
	else:
		theta[CS-1]=float(var)
	replot(params_float_zoomed[0],params_float_zoomed[1],N)
	
theta_slider= ScaleNew(root, from_=theta_min, to=theta_max, orient=HORIZONTAL, background='white',highlightbackground='white',  sliderlength=5, width=10, length=258, resolution=abs(theta_max-theta_min)/360, sliderrelief=FLAT, troughcolor="#ffddcc", command=update_theta)
theta_slider.grid(row =5, column=41, columnspan=20, padx=(0,30), sticky=W)

def update_kl(var):
	global kl
	kl=float(var) 
	replot(params_float_zoomed[0],params_float_zoomed[1],N)
	
kl_slider= ScaleNew(root, from_=kl_min, to=kl_max, orient=HORIZONTAL, background='white',highlightbackground='white',  sliderlength=5, width=10, length=258, resolution=abs(kl_max-kl_min)/200, sliderrelief=FLAT, troughcolor="#ffddcc", command=update_kl)
kl_slider.grid(row =9, column=41, columnspan=20, padx=(0,30), sticky=W)
kl_slider.set(kl)

def update_nu(var):
	global nu
	global wp
	global cc
	global nu_sliderVar
	
	nu_sliderVar[CS-1]=var
	
	nu[CS-1]=wp[0]*10**float(var) 
	
	nu_value_label2 = Text(root, width=30, height=2, borderwidth=0,background=root.cget("background"))
	nu_value_label2.tag_configure("subscript", offset=-4)
	nu_value_label2.insert("insert", u"\u03bd: 10^"+ '{0:.2f}'.format(float(var))+" *"+u"\u03c9", "", "p1", "subscript")
	nu_value_label2.configure(state="disabled")
	nu_value_label2.grid(row=7, column=62,columnspan=12,padx=(0,0), sticky=W)
	replot(params_float_zoomed[0],params_float_zoomed[1],N)
nu_slider= ScaleNew(root, from_=nu_min_exponent, to=nu_max_exponent+4, orient=HORIZONTAL, background='white',highlightbackground='white',  sliderlength=5, width=10, length=272, resolution=abs(nu_max_exponent-nu_min_exponent)/nu_res, sliderrelief=FLAT, troughcolor="#ffddcc", command=update_nu)
nu_slider.grid(row =5, column=62, columnspan=12, padx=(0,10), sticky=W)
nu_slider.set(nu_sliderVar[CS-1])


def update_wp(var):
	global wc
	global wp
	wp[CS-1]=float(var) 
	
	replot(params_float_zoomed[0],params_float_zoomed[1],N)
wp_slider= ScaleNew(root, from_=wp_min, to=wp_max, orient=HORIZONTAL, background='white',highlightbackground='white',  sliderlength=5, width=10, length=272, resolution=abs(wp_max-wp_min)/200, sliderrelief=FLAT, troughcolor="#ffddcc", command=update_wp)
wp_slider.grid(row =9, column=62, columnspan=12, padx=(0,10), sticky=W)	
wp_slider.set(wp[CS-1])

def update_sliders():
	global wp_slider
	global theta_slider
	global nu_slider
	global nu_sliderVar
	
	nu_slider.set(nu_sliderVar[CS-1])
	wp_slider.set(wp[ CS-1])
	theta_slider.set(theta[CS-1])




ax_zoom = fig.add_axes([0.12, 0.1, 0.78, 0.03]) #position and dimensions for slider
zoom_min=-2
zoom_max=2
zoom = Slider(ax_zoom, 'Zoom', zoom_min, zoom_max, valinit=0)

zoom_coord_entry=Entry(root, width=7, relief=FLAT, highlightbackground="blue", highlightthickness=2)
zoom_coord_entry.grid(row=39, column=2,columnspan=1,padx=(0,10), sticky=W)
zoom_coord_label= Label(root, width=14,text='Zoom to X,Y:', background='white', font=("Courier",10))
zoom_coord_label.grid(row=39, column=0,columnspan=2,padx=(16,0), sticky=W)

scroll_check=0

def update_zoom(val):
	global params_float_zoomed
	global wmin
	global wmax
	global scroll_check
	if(scroll_check==0):
		zcs=zoom_coord_entry.get() #Zoom Coordinate String = zcs
		rr = re.findall("[+-]?\d*[\.]?\d*(?:(?:[eE])[+-]?\d+)?", zcs)
		i=0
		x_center=(params_float_zoomed[1]+params_float_zoomed[0])/2
		y_center=(params_float_zoomed[3]+params_float_zoomed[2])/2
		zc_float=[x_center,y_center]
		for numbers in rr:
			if(numbers!="" and i<2):
				zc_float[i]=float(numbers)
				i+=1
		x_coord=zc_float[0]
		y_coord=zc_float[1]
		#print([x_center, y_center, x_coord,y_coord])
		#originally i wanted to test for ValueErrors, BUT when there is a wrong input in the entry box, then the for-loop does not go into the if-condition, and will not replace the preset (zoom to center coordinates)
		
		
		z = zoom.val #reminder that this is some matplotlib thing, so it doesnt work like the tkinter slider! .val is giving an error on tkinter sliders!!
		#zoom to entried coordinates 
		
		params_float_zoomed=[x_coord-(x_coord-params_float[0])/10**z,x_coord+(params_float[1]-x_coord)/10**z,y_coord-(y_coord-params_float[2])/10**z,y_coord+(params_float[3]-y_coord)/10**z]
		ax.axis(params_float_zoomed)
		wmin=params_float_zoomed[0]
		wmax=params_float_zoomed[1]
		replot(params_float_zoomed[0],params_float_zoomed[1],N)
	
zoom.on_changed(update_zoom) #when value is changed, then on_changed() is called


def update_zoom_mouse(event):
	#mouse has 3 events: delta, x, y 
	#delta returns a constant value that indicates the scrolling of the mouswheel. every "step" one gets the same value, with the sign depending on the scroll direction. 
	#when i wrote this, the constant value for the mouse in use was 120, or -120
	#left border of graph:x=80, right border: x=576 (do not change size of graph!!)
	#upper border of graph:y=10, lower border: y=360
	global params_float_zoomed
	global wmin
	global wmax
	global pos_coord
	global scroll_check
	r_border=576 #right
	l_border=80 #left
	t_border=10 #top
	b_border=360 #bottom
	x_pos= event.x
	y_pos=event.y
	delta=event.delta
	#mouse position in the units of the graph
	x_coord = params_float[0]+(params_float[1]-params_float[0])*(x_pos-l_border)/(r_border-l_border)
	y_coord = params_float[2]+(params_float[3]-params_float[2])*(y_pos-b_border)/(t_border-b_border)
	x_coord_zoomed = params_float_zoomed[0]+(params_float_zoomed[1]-params_float_zoomed[0])*(x_pos-l_border)/(r_border-l_border)
	y_coord_zoomed = params_float_zoomed[2]+(params_float_zoomed[3]-params_float_zoomed[2])*(y_pos-b_border)/(t_border-b_border)
	
	z = zoom.val
	z=z +(zoom_max-zoom_min)/12000*delta #increase z by a (delta/....)-th of the scale range
	params_float_zoomed=[x_coord_zoomed-(x_coord-params_float[0])/10**z,x_coord_zoomed+(params_float[1]-x_coord)/10**z,y_coord_zoomed-(y_coord-params_float[2])/10**z,y_coord_zoomed+(params_float[3]-y_coord)/10**z]
	#zoomed coordinates as reference for position, original coordinates as reference for zoom (they are at all times multiplied by the zoom factor! putting the zoomed coordinates in there will zoom 'twice')
	#ax.axis([x_coord-(x_coord-params_float[0])/10**z,x_coord+(params_float[1]-x_coord)/10**z,y_coord-(y_coord-params_float[2])/10**z,y_coord+(params_float[3]-y_coord)/10**z])
	ax.axis(params_float_zoomed)
	wmin=params_float_zoomed[0]
	wmax=params_float_zoomed[1]
	replot(params_float_zoomed[0],params_float_zoomed[1],N)
	#only when the zoom scale is manually operated, it should zoom. when the wheel is scrolled, it should update only the position of the scale.
	#but updating the position of the scale, changes its value, thus zoom.on_changed() is called (see above the update_zoom_mouse function)
	#in order to have only the scale value updated, and not the function executed, i introduced a "scroll check"
	scroll_check=1
	zoom.set_val(z) #has to be after zoomed parameters!!
	scroll_check=0
canvas.get_tk_widget().bind("<MouseWheel>", update_zoom_mouse)


def update_param():
	
	parameters_entry=parameters.get()
	#rr = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", parameters_entry)
	rr = re.findall("[+-]?\d*[\.]?\d*(?:(?:[eE])[+-]?\d+)?", parameters_entry)
	global params_float 
	global params_float_zoomed
	
	i=0
	for numbers in rr:
		if(numbers!="" and i<4):
			params_float[i]=float(numbers)
			i+=1
	params_float_zoomed = params_float
	global zoom
	zoom.reset() #puts the zoom scale back to 1
	
	ax.axis(params_float)
	replot(params_float[0],params_float[1],N) #needed so the length of the plotted arrays are adjusted
	fig.canvas.draw_idle() #needed additionally to update the scale *even if* nothing is plotted!
	
parameters.bind("<Return>", (lambda event: update_param()))


def Plot_CP_ellipse(event):
			global wmin
			global wmax
			global params_float
			global params_float_zoomed
			global getcoord
			global lineCP1
			global lineCP2
			global plotCPcondition
			global CPButtonCheck
			global CPsave 
			global phi
			global N 
			global fig
			global canvas 
			global wc 
			global kl
			global wp
			
			
			r_border=576 #right
			l_border=80 #left
			t_border=10 #top
			b_border=360 #bottom
			x_pos= event.x
			y_pos=event.y
			delta=event.delta
			#mouse position in the units of the graph
			
			if(CPButtonCheck==1):	#avoid multiple calls
				x_coord_CP = params_float_zoomed[0]+(params_float_zoomed[1]-params_float_zoomed[0])*(x_pos-l_border)/(r_border-l_border)
				y_coord_CP = params_float_zoomed[2]+(params_float_zoomed[3]-params_float_zoomed[2])*(y_pos-b_border)/(t_border-b_border)
				#wp=1
				w=x_coord_CP #not globally defined, so i just made it a single valued variable
				swapper=(wp[0]-w)/abs(wp[0]-w)
				
				
				
				k1=kl
				k2=kl
				k0=kl
				kxx=kl
				kxy=0
				kzz=kl
				
				cc=0
				
				
				while cc<=NCarriers-1:
					mrel=wp[0]**2/wp[cc]**2 #=mc/m1**2 since wp ~ 1/sqrt(m)
					swapper=(wp[cc]-w)/abs(wp[cc]-w) #this allows to swap between k1 and k2 definitions at wp, thus giving R+ and R- rather than R1 and R2, in regards to Paliks definitions 
					#see page 1230 in palik/furdyna:
					k1+=-kl*(2*wp[cc]**2*(w*(w+1.j*nu[cc])-wp[cc]**2))/(w*(2*(w+1.j*nu[cc])*(w*(w+1.j*nu[cc])-wp[cc]**2)-w*(wc*mrel)**2*np.sin(theta[cc]/360*2*np.pi)**2+swapper*(wc*mrel)*(w**2*(wc*mrel)**2*np.sin(theta[cc]/360*2*np.pi)**4+4*(w*(w+1.j*nu[cc])-wp[cc]**2)**2*np.cos(theta[cc]/360*2*np.pi)**2)**0.5))
					k2+=-kl*(2*wp[cc]**2*(w*(w+1.j*nu[cc])-wp[cc]**2))/(w*(2*(w+1.j*nu[cc])*(w*(w+1.j*nu[cc])-wp[cc]**2)-w*(wc*mrel)**2*np.sin(theta[cc]/360*2*np.pi)**2-swapper*(wc*mrel)*(w**2*(wc*mrel)**2*np.sin(theta[cc]/360*2*np.pi)**4+4*(w*(w+1.j*nu[cc])-wp[cc]**2)**2*np.cos(theta[cc]/360*2*np.pi)**2)**0.5))
					#see page 1205 in palik/furdyna::
					kxx+=(1.0j/w)*wp[cc]**2*kl*((nu[cc]-1.0j*w)/((nu[cc]-1.0j*w)**2+(wc*mrel)**2))
					kxy+=(1.0j/w)*wp[cc]**2*kl*((wc*mrel)/((nu[cc]-1.0j*w)**2+(wc*mrel)**2))
					kzz+=(1.0j/w)*wp[cc]**2*kl*(1.0/(nu[cc]-1.0j*w))
					cc+=1
				cc=0

				try:
					ExOverEy1_CP=(kzz/kxy)*(kxx-k1)*np.cos(theta[cc]/360*2*np.pi)/(kzz - k1*np.sin(theta[cc]/360*2*np.pi)**2)
					ExOverEy2_CP=(kzz/kxy)*(kxx-k2)*np.cos(theta[cc]/360*2*np.pi)/(kzz - k2*np.sin(theta[cc]/360*2*np.pi)**2)
				except ZeroDivisionError:
					ExOverEy1_CP=0
					ExOverEy2_CP=0
				
				
				#could not manage to overlay 2 independent subplots. if they are both defined on (111), then the second one will be ignored, or rather made equal to the first...so i could not introduce separate axes
				#since i have to work in the same coordinates as the other plots, i have to apply weird scalings to the ellipses to make them look correct. if not scaled correctly, a circle will not look like a circle in these coordinates!
				getcoord_CP=[x_coord_CP,y_coord_CP]
				w_circ=1/(r_border-l_border)*(params_float_zoomed[1]-params_float_zoomed[0])
				h_circ=1/(t_border-b_border)*(params_float_zoomed[3]-params_float_zoomed[2])
				if(abs(ExOverEy1_CP.imag)<1):
					x_ellipse1=ExOverEy1_CP.imag*100*w_circ*np.cos(phi)+x_coord_CP
					y_ellipse1=100*h_circ*np.sin(phi)+params_float_zoomed[2]+0.5*(params_float_zoomed[3]-params_float_zoomed[2])*(0-b_border)/(t_border-b_border)
				else:
					x_ellipse1=100*w_circ*np.cos(phi)+x_coord_CP
					y_ellipse1=1/ExOverEy1_CP.imag*100*h_circ*np.sin(phi)+params_float_zoomed[2]+0.5*(params_float_zoomed[3]-params_float_zoomed[2])*(0-b_border)/(t_border-b_border)
				
				if(abs(ExOverEy2_CP.imag)<1):
					x_ellipse2=ExOverEy2_CP.imag*100*w_circ*np.cos(phi)+x_coord_CP
					y_ellipse2=100*h_circ*np.sin(phi)+params_float_zoomed[2]+0.5*(params_float_zoomed[3]-params_float_zoomed[2])*(0-b_border)/(t_border-b_border)
				else:
					x_ellipse2=100*w_circ*np.cos(phi)+x_coord_CP
					y_ellipse2=1/ExOverEy2_CP.imag*100*h_circ*np.sin(phi)+params_float_zoomed[2]+0.5*(params_float_zoomed[3]-params_float_zoomed[2])*(0-b_border)/(t_border-b_border)
				
				
				if(plotCPcondition==0):
					lineCP1, = ax.plot(x_ellipse1, y_ellipse1) 
					lineCP2, = ax.plot(x_ellipse2, y_ellipse2) 
					plotCPcondition=1
				if(plotCPcondition>0):
					lineCP1.set_ydata(y_ellipse1)
					lineCP2.set_ydata(y_ellipse2)
					lineCP1.set_xdata(x_ellipse1)
					lineCP2.set_xdata(x_ellipse2)
				
				fig.canvas.draw_idle()
		

mainloop()