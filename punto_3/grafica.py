import matplotlib.pyplot as plt
import numpy as np
import corner

data = np.loadtxt('data.txt')
best_alpha = np.median(data[:, 0])
best_betta = np.median(data[:, 1])
best_gamma = np.median(data[:, 2])
best_delta = np.median(data[:, 3])
print "El valor mas probable de alpha es ", best_alpha
print "El valor mas probable de beta es ", best_betta
print "El valor mas probable de gamma es ", best_gamma
print "El valor mas probable de delta es ", best_delta
# No pude poner la media y los errores porque al poner "show_titles=True" me mandaba un error que no supe resolver
# Si su corazon es bondadoso, podria intentar meterlo para ver si es solo mi computador? Gracias!
figure = corner.corner(data, labels=[r"$\alpha$", r"$\beta$", r"$\gamma$", r"$\delta$"])
plt.savefig('grafica.pdf')
