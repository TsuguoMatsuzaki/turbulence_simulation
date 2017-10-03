FC = mpifrtpx
TARGET = a.out

SRCDIR = .
OBJDIR = ./.obj
MODDIR = ./.mod
LSTDIR = ./.lst

FFLAGS = -M$(MODDIR) -Qo$(LSTDIR)/$(@F:.o=.lst),t -Kfast,openmp,preex
LIBS = -lfftw3 -lm
LDFLAGS = -L/home/apps/fftw/3.3.3/lib64

OBJS =\
	$(OBJDIR)/common.o\
	$(OBJDIR)/exchange.o\
	$(OBJDIR)/fft.o\
	$(OBJDIR)/forcing.o\
	$(OBJDIR)/integration.o\
	$(OBJDIR)/io_fields.o\
	$(OBJDIR)/main.o\
	$(OBJDIR)/mpi_common.o\
	$(OBJDIR)/nonlinear.o\
	$(OBJDIR)/output.o\
	$(OBJDIR)/param.o\
	$(OBJDIR)/spectrum.o\
	$(OBJDIR)/timer.o

all : init $(TARGET)

$(TARGET) : $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -o $@ -c $<

$(OBJDIR)/common.o :\
	$(OBJDIR)/param.o\
	$(OBJDIR)/mpi_common.o

$(OBJDIR)/exchange.o :\
	$(OBJDIR)/param.o\
	$(OBJDIR)/mpi_common.o

$(OBJDIR)/fft.o :\
	$(OBJDIR)/param.o\
	$(OBJDIR)/exchange.o\
	$(OBJDIR)/timer.o

$(OBJDIR)/forcing.o :\
	$(OBJDIR)/param.o\
	$(OBJDIR)/mpi_common.o\
	$(OBJDIR)/common.o

$(OBJDIR)/integration.o :\
	$(OBJDIR)/param.o\
	$(OBJDIR)/mpi_common.o\
	$(OBJDIR)/common.o\
	$(OBJDIR)/forcing.o\
	$(OBJDIR)/nonlinear.o\
	$(OBJDIR)/timer.o

$(OBJDIR)/io_fields.o :\
	$(OBJDIR)/param.o\
	$(OBJDIR)/mpi_common.o

$(OBJDIR)/main.o :\
	$(OBJDIR)/param.o\
	$(OBJDIR)/mpi_common.o\
	$(OBJDIR)/common.o\
	$(OBJDIR)/io_fields.o\
	$(OBJDIR)/output.o\
	$(OBJDIR)/fft.o\
	$(OBJDIR)/spectrum.o\
	$(OBJDIR)/nonlinear.o\
	$(OBJDIR)/integration.o\
	$(OBJDIR)/timer.o

$(OBJDIR)/mpi_common.o :\
	$(OBJDIR)/param.o

$(OBJDIR)/nonlinear.o :\
	$(OBJDIR)/param.o\
	$(OBJDIR)/mpi_common.o\
	$(OBJDIR)/common.o\
	$(OBJDIR)/fft.o\
	$(OBJDIR)/timer.o

$(OBJDIR)/output.o :\
	$(OBJDIR)/param.o\
	$(OBJDIR)/mpi_common.o\
	$(OBJDIR)/common.o

$(OBJDIR)/spectrum.o :\
	$(OBJDIR)/param.o\
	$(OBJDIR)/mpi_common.o\
	$(OBJDIR)/common.o

$(OBJDIR)/timer.o :\
	$(OBJDIR)/mpi_common.o

init :
	@if [ ! -d $(OBJDIR) ]; then\
		mkdir $(OBJDIR);\
	fi
	@if [ ! -d $(MODDIR) ]; then\
	   mkdir $(MODDIR);\
	fi
	@if [ ! -d $(LSTDIR) ]; then\
	   mkdir $(LSTDIR);\
	fi
clean :
	@rm -rf $(TARGET) $(OBJDIR)/*.o $(MODDIR)/*.mod $(LSTDIR)/*.lst

cleanall :
	@rm -rf $(TARGET) $(OBJDIR) $(MODDIR) $(LSTDIR)
