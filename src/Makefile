LIBNAME=growth
LD=ld

OBJS=	carmasub.o cgamma.o elliptic.o

$(LIBNAME).$(SHLIBEXT):	$(OBJS)
	$(LD) $(SHLIBLDFLAGS) -o $@ $(OBJS)

clean:
	@rm -f *.o *.$(SHLIBEXT)

realclean: clean
