CC = g++
RVERSION=2011_3
CPPFLAGS = -Wall -g -O3 -DNDEBUG
NOFLG = -DNDEBUG

INDEXOBJS = Tools.o BitRank.o HuffWT.o TBWTIndex.o

all: builder test traverse checkresult queries

traverse: traverse.o $(INDEXOBJS)
	$(CC) $(CPPFLAGS) -o traverse traverse.o $(INDEXOBJS)

queries: queries.o $(INDEXOBJS)
	$(CC) $(CPPFLAGS) -o queries queries.o $(INDEXOBJS)

test: test.o $(INDEXOBJS)
	$(CC) $(CPPFLAGS) -o test test.o $(INDEXOBJS)

checkresult: checkresult.o $(INDEXOBJS)
	$(CC) $(CPPFLAGS) -o checkresult checkresult.o $(INDEXOBJS)

builder: builder.o TBWTBuilder.o $(INDEXOBJS)
	$(CC) $(CPPFLAGS) -o builder builder.o TBWTBuilder.o $(INDEXOBJS)

depend:
	g++ -MM *.cpp > dependencies.mk

new_release:
	@echo "Creating new release, version $(RAVERSION), into subdir tbwt_$(RAVERSION)"
	@echo "Press <enter> to continue or <ctrl-c> to abort:"
	@read unused
	@make -C $(LIBCDSPATH) clean
	@make -C $(LIBRLCSAPATH) clean
	@mkdir tbwt_$(RAVERSION)
	@cp Makefile tbwt_$(RAVERSION)/
	@cp *.cpp *.h tbwt_$(RAVERSION)/
	@make -C tbwt_$(RAVERSION) depend
	@tar zcf tbwt_$(RAVERSION).tar.gz --exclude-vcs tbwt_$(RAVERSION)
	@rm -rf tbwt_$(RAVERSION)
	@echo "Done, new file tbwt_$(RAVERSION).tar.gz created."

clean:
	rm -f core *.o *~ builder traverse test checkresult queries tconvert resultconvert

shallow_clean:
	rm -f core *.o *~ builder traverse test

include dependencies.mk
