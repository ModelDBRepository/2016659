
CC= g++ 
#CPPFLAGS=-g  -Wall -std=gnu++1z
CPPFLAGS= -O3 -std=gnu++1z

#INCS=...

LIBSRCS= util.cpp simIO.cpp

LIBOBJS= $(LIBSRCS:.cpp=.o) 

all: simulator 

simulator:  $(LIBOBJS) simulator.cpp  
	$(CC) $(CPPFLAGS) simulator.cpp -o simulator $(LIBOBJS) $(LIBS)

clean: 
	rm *.o *.d  *.a 

#%.o: %.cpp
	#$(CC)  -MD $(CPPFLAGS) -c -o $@ $<

%.o: %.cpp 
	$(CC)  $(CPPFLAGS) -c -o $@ $<

