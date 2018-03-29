GraphPlanar.so: Makefile

PLANARITYSRC=graphs/graphIsolator.c graphs/graphStructure.c graphs/graphEmbed.c graphs/graphNonplanar.c graphs/graphTests.c   graphs/stack.c graphs/graphIO.c graphs/graphPreprocess.c graphs/listcoll.c

PLANARITYOBJ=$(PLANARITYSRC:.c=.o)

PYTHONPREFIX=$(shell python -c "import sys; print sys.prefix")
PYTHONVER=$(shell python -c "import sys; print '%i.%i'%(sys.version_info[:2])")
PYTHONINCLUDE=$(PYTHONPREFIX)/include/python$(PYTHONVER)

GraphPlanar.so: GraphPlanar.o $(PLANARITYOBJ)
	gcc -shared GraphPlanar.o $(PLANARITYOBJ) -fPIC -o $@

%.o : %.c
	gcc -c -I $(PYTHONINCLUDE) -fPIC -I graphs/ $< -o $@


clean:
	rm -fv $(PLANARITYOBJ) GraphPlanar.o GraphPlanar.so
