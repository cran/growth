LIBNAME=growth

OBJS=	carmasub.o cgamma.o elliptic.o resid.o

$(LIBNAME).$(SHLIBEXT):	$(OBJS)
	$(SHLIBLD) $(SHLIBLDFLAGS) $(FLIBS) -o $@ $(OBJS)

clean:
	@rm -f *.o *.$(SHLIBEXT)

realclean: clean
