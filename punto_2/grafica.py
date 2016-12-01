import matplotlib.pyplot as plt
import numpy as np
import corner

data = np.loadtxt('data.txt')
best_alpha = np.median(data[:, 0])
best_log_Ms = np.median(data[:, 1])
print "El valor mas probable de alpha es ", best_alpha, ", y el valor mas probable del logaritmo base 10 de la masa del sol es ", best_log_Ms
# No pude poner la media y los errores porque al poner "show_titles=True" me mandaba un error que no supe resolver
# Si su corazon es bondadoso, pordria intentar meterlo para ver si es solo mi computador? Gracias!
figure = corner.corner(data, labels=[r"$\alpha$", r"$\log M_s$"])
plt.savefig('grafica.pdf')
