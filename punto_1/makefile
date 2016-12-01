grafica.pdf : data.txt grafica.py
	python grafica.py

data.txt : salida.x
	./salida.x > data.txt

salida.x : sismo.c
	gcc sismo.c -o salida.x -lm

clean:
	rm -f salida.x data.txt
