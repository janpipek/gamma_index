libgamma.so: gamma.c
	gcc -O3 gamma.c -lm -fPIC -shared -o libgamma.so

test.exe: gamma.c
	gcc -O3 -o test.exe -lm gamma.c
	./test.exe

clean:
	rm libgamma.so test.exe