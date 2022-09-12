APP				= main
CC 				= clang++
CFLAGS			= -lm -std=c++0x	# -g For debugging.
PYTHON_INCLUDE 	= /Library/Developer/CommandLineTools/SDKs/MacOSX11.1.sdk/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7

FILENAME1		= ./transforms

default: compile
	
compile:
	$(CC) $(CFLAGS) -I kissfft kissfft/kiss_fft.c -c $(FILENAME1).cpp
	$(CC) $(CFLAGS) -lsndfile -I $(PYTHON_INCLUDE) -lpython2.7 $(APP).cpp -o $(APP) $(FILENAME1).o kiss_fft.o
	@echo "Build done."
	
run:
	./main samples/v2-please.wav classifier database.txt

clean:
	rm $(FILENAME1).o
	rm $(APP)
