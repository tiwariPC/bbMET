import matplotlib.pyplot as plt

mphi=[50.,100.,350.,400.,500.,1000.]
sr1eff=[0.02110994, 0.03111132, 0.06311873, 0.0679374, 0.07349698, 0.09729622]
sr2eff=[0.00339505, 0.00392721, 0.00878031, 0.0144919, 0.01378335, 0.02291534]

plt.plot(mphi,sr1eff,'o-',label="SR1")
plt.plot(mphi,sr2eff,'o-',label="SR2")
plt.xlabel(r'$M_{\phi}$')
plt.ylabel("Signal Efficiency (fraction)")
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlim(40,1200)
plt.title("pseudoscalar, $M_{\chi}$ = 1 GeV")
ticks=[50.,100.,350.,500.,1000.]
plt.xticks(ticks,ticks)
plt.grid()
#plt.show()
plt.savefig("Pseudo.png")
