HOST	      = DESKTOP-FK

DEST	      = /home/koikef/bin      

EXTHDRS	      =

HDRS	      =

FC	      = gfortran

FFLAGS	      = -g

LDFLAGS	      = -g

LIBS	      =

LINKER	      = gfortran

MAKEFILE      = Makefile

OBJS	      = Main.o \
		VPA-Lib.o \
		RICBES.o \
		Potential.o \
		DLSODALIB.o \
		WKB.o \
		cgcoef.o \
		legendre.o

PRINT	      =    

PRINTER	      =      

PROGRAM	      = difelsccs

SRCS	      = Main.f90 \
		VPA-Lib.f90 \
		RICBES.f \
		Potential.f90 \
		DLSODALIB.f \
		WKB.f90 \
		cgcoef.f \
		legendre.f90

PORTABLE      = VPA.FORT

all:		$(PROGRAM)

$(PROGRAM):	$(SRCS)
		@$(FC) $(FFLAGS) -o $(PROGRAM) $(SRCS)
		@echo -n "Loading $(PROGRAM) ... "
		@echo "done"

clean:;		@rm -f $(OBJS)

depend:;	@mkmf -f $(MAKEFILE) PROGRAM=$(PROGRAM) DEST=$(DEST)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@cp $(PROGRAM).$(HOST) $(DEST)/$(PROGRAM)
		 
print:;		@$(PRINT) $(HDRS) $(SRCS)

list:		$(SRCS)
		@$(PRINT) $?
		@echo '' >> list
		@echo Listed files are: >> list
		@echo $? >> list
		@echo `date` >> list

program:        $(PROGRAM)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

combine:;       @cat $(SRCS)|tr a-z A-Z > $(PORTABLE)

$(DEST)/$(PROGRAM): $(SRCS) $(LIBS) $(HDRS) $(EXTHDRS)
		@make -f $(MAKEFILE) DEST=$(DEST) install

RCS:		
		@mkdir ./RCS

ci:;		@make RCS
		@ci -l $(SRCS)

co:;		@co -q -l $(SRCS)

register:;	@make RCS
		@ci $(SRCS)
		@co -q -l $(SRCS)

$(SRCS):;	@co $(SRCS)
