#Object files

OBJ = main.o OptionParser.o properties.o TableCompiler.o utilities.o map.o Output.o ReadData.o prediction.cpp SNP.o

#Header files

HDR = main.h OptionParser.h properties.h TableCompiler.h utilities.h common.h map.h Output.h ReadData.h prediction.h SNP.h


default: help
install: prediction

CPP= g++ -c -o $@ $<
LINK=g++ -o $@

#All the rules

%.o : %.cpp ${HDR}
	${CPP}

#The final step
prediction : ${OBJ}
	${LINK} ${OBJ}

clean:
	/bin/rm -f *.o prediction
	/bin/rm -f *.*~

help:
	@echo "make clean    => Clean compiled objs and executables"
	@echo "make install  => Compile the program"
