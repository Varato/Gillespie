import numpy as np
import matplotlib.pyplot as plt

# reaction rate constants
k1 = 2.
k2 = 1.
k3 = 1.

# the system is described as a dictionary
system = dict()

state = 0 # global state
t = 0     # global time
state_record = []
t_record = []
# trans_record = []
T0 = 0
T1 = 0

def init(Cs):
	global state, t, state_record, t_record, system, T0, T1
	state = 0
	t = 0
	state_record = [0]
	t_record = [0]
	system = {0:[k1*Cs], 1:[k2, k3]}
	T0=0
	T1=1

def one_MC_step():
	global state, t, state_record, t_record, trans_record, T0, T1
	paths = system[state]

	Rtot = sum(paths)

	# sampling 1 transition
	trans_time = [-np.log(np.random.random())/k for k in paths]
	dt = min(trans_time)
	
	if state==0:
		T0 = T0+dt
	elif state==1:
		T1 = T1+dt
	state = {0:1, 1:0}[state]
	t = t + dt
	t_record.append(t)
	state_record.append(state)

def main(MC_steps, Cs):
	init(Cs)
	for i in range(MC_steps):
		one_MC_step()

	T = t_record[-1]
	p0 = T0/T
	p1 = T1/T
	v = p1*k3
	return (p0, p1, v)

def paint_state(ax, Cs):
	# draw the signal
	signal=[ [t_record[0], state_record[0]] ]
	for j in range(1,len(state_record)):
		signal.append([t_record[j], state_record[j-1]])
		signal.append([t_record[j], state_record[j]])

	signal = np.array(signal)
	# for i in range(0, len(signal)-1, 2):
	# 	plt.plot( [signal[i][0], signal[i+1][0]], [signal[i][1], signal[i+1][1]], 'r-' )

	ax.plot(signal[:,0], signal[:,1], '-', color="blue")
	ax.axis([0, 10, -0.1, 1.1])
	ax.set_yticks([0, 1])
	# ax.set_ylabel("state", fontsize=30)
	# ax.set_xlabel("time", fontsize=30)
	ax.tick_params(labelsize=20)
	ax.set_title("Cs={}".format(Cs),fontsize=25)

def Analytic_curve(ax):
	Cs = np.linspace(0, 50, 500)
	KH=(k2+k3)/k1

	def MM(Cs):
		return k3*Cs/(Cs+KH)

	ax.plot(Cs, MM(Cs), label="analytic curve")


if __name__=="__main__":
	fig1, ax1 = plt.subplots(figsize=[10, 10])
	fig2, ((ax21,ax22),(ax23,ax24))=plt.subplots(2,2, figsize=[10,10],sharex=True, sharey=True)
	fig2.text(0.5,0.03,"time",ha="center", fontsize=25)
	fig2.text(0.03,0.5,"states",va="center", rotation="vertical",fontsize=25)
	fig2.text(0.5,0.96,"state-time dependence",ha="center", fontsize=25)
	Cs_vec = [0.1, 0.5, 1, 2, 5, 10, 20, 30, 50]
	v_ave = []
	v_std = []
	cycles = 10
	for Cs in Cs_vec:
		print("Current Cs = ", Cs)
		v_record = []
		for c in range(cycles):
			(p0, p1, v) = main(10000, Cs)
			v_record.append(v)
		v_record = np.array(v_record)
		v_ave.append(np.mean(v_record))
		v_std.append(np.std(v_record))
		if Cs==0.5:
			paint_state(ax21, Cs)
		if Cs==2:
			paint_state(ax22, Cs)
		if Cs==10:
			paint_state(ax23, Cs)
		if Cs==50:
			paint_state(ax24, Cs)
		print(p0, p1)
		print( (k2+k3)/(k1*Cs+k2+k3) )
	print(v_std)
	Analytic_curve(ax1)
	ax1.errorbar(Cs_vec, v_ave, yerr=v_std, color="red", fmt="o", label="Gillespie result")
	ax1.set_ylim([0,k3+k3*0.1])
	ax1.set_xlabel("Cs(a.u.)", fontsize=25)
	ax1.set_ylabel("$v$(a.u.)", fontsize=25)
	ax1.legend(loc="lower right", fontsize=25)
	ax1.tick_params(labelsize=20)
	ax1.set_title("M.M. curve @ k1={}, k2={}, k3={}".format(k1,k2,k3),fontsize=25)

	fig1.savefig("/Users/xinchen/Dropbox/Biophysics_final_project/report/img/MM2.png", dpi=300)
	fig2.savefig("/Users/xinchen/Dropbox/Biophysics_final_project/report/img/state2.png", dpi=300)
	# plt.show()




