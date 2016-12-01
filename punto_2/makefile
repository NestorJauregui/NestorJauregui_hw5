grafica.pdf : data.txt grafica.py
	python grafica.py

data.txt : salida.x
	./salida.x > data.txt

salida.x : sol.c
	gcc sol.c -o salida.x -lm

clean:
	rm -f salida.x data.txt
