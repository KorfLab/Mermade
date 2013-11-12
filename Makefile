############################
# Makefile for motifamatic #
############################

LIB = -lm

APP = motifamatic
SRC = motifamatic.c
OBJ = motifamatic.o

DATE = $(shell date +\%Y-\%m-\%d)

###########
# Targets #
###########

default:
	make gcc

$(APP): $(OBJ)
	$(CC) -o $(APP) $(CFLAGS) $(OBJ) $(LIB)

clean:
	rm -f *.o $(APP)
	
tar:
	mkdir mermade-1.0
	cp README Makefile motifamatic.c sequence_converter.pl de-barcode.pl kmer_counter.pl kmer_selector.pl mermade.pl motif_expander.pl db_creator.pl run_mermade.pl reporter.pl README.pdf mermade-1.0
	tar -zvcf mermade-1.0.tar.gz mermade-1.0
	rm -rf mermade-1.0

#################
# Architectures #
#################

gcc:
	make $(APP) CC="gcc" CFLAGS="-O2 -Wall -Werror"

###################
# Inference Rules #
###################

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<


