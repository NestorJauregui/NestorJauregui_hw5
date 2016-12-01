import matplotlib.pyplot as plt
import numpy as np
import corner

data = np.loadtxt('data.txt')
best_x = np.median(data[:, 0])
best_y = np.median(data[:, 1])
print "El valor mas probable de las coordenadas del epicentro son (", best_x, ",", best_y, ")."
# No pude poner la media y los errores porque al poner "show_titles=True" me mandaba un error que no supe resolver
# Si su corazon es bondadoso, podria intentar meterlo para ver si es solo mi computador? Gracias!
figure = corner.corner(data, labels=[r"$x$", r"$y$"])
plt.savefig('grafica.pdf')
