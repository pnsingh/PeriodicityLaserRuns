CPP = g++
CC = gcc
CFLAGS = $(shell root-config --cflags) 
LIBS = $(shell root-config --glibs) -lMinuit -lMinuit2 -lFoam -lfftw3 
GLIBS = 
GLIBS += 
OBJECTS = TestFft.o
HEADERS = 

ALL : TestFft.exe
	echo "Listo!"

TestFft.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o TestFft.exe $(LIBS) $(GLIBS) $(CFLAGS)

TestFft.o : TestFft.cc $(HEADERS)
	$(CPP) -c TestFft.cc -o TestFft.o $(CFLAGS) $(LIBS)

clean:
	rm -f *~ *.o *.exe 
