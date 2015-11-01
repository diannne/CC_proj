IDIR =../include
CC=g++
CFLAGS=-I$(IDIR)

ODIR=obj
LDIR =../lib

LIBS=-lfreeimage

_DEPS = ColorTransfer.h KImage.h Histogram.h MatrixModule.h Color2Gray.h Color_Models_Conversions.h	PdfTransfer.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = Color2Gray.o ColorTransfer.o KImage.o Histogram.o Main.o MatrixModule.o  PdfTransfer.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

NPdfTransfer: $(OBJ)
	g++ -L$(LDIR) -o $@ $^ $(CFLAGS) $(LIBS) 

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~.
