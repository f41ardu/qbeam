all: qbeam 

f77FLAGS = -r8R

OBJ = angles.o blc8.o blr8.o brtheta.o c.o calc.o cross.o dot.o ep-src-s.o lin8c.o lin8r.o myeb.o myexp.o qbeam4.o svbksb.o tensor.o gikn.o

LFLAGS = -lm -lc 

qbeam: $(OBJ)
	gfortran -o qbeam $(OBJ) $(LFLAGS)

clean:
	rm -f qbeam $(OBJ)

.f.o:
	gfortran $(f77FLAGS) -c $<

FRC:
	.SUFFFIXES: .f

# DEPENDENCIES

