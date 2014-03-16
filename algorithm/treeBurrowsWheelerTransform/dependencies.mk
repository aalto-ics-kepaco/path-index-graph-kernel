BitRank.o: BitRank.cpp BitRank.h Tools.h
HuffWT.o: HuffWT.cpp HuffWT.h BitRank.h Tools.h
TBWTBuilder.o: TBWTBuilder.cpp TBWTBuilder.h SimpleForest.h BitRank.h \
  Tools.h BlockArray.h
TBWTIndex.o: TBWTIndex.cpp TBWTIndex.h BitRank.h Tools.h HuffWT.h \
  BlockArray.h
Tools.o: Tools.cpp Tools.h
builder.o: builder.cpp SimpleForest.h BitRank.h Tools.h TBWTBuilder.h \
  BlockArray.h HuffWT.h
checkresult.o: checkresult.cpp NaiveForest.h Tools.h
queries.o: queries.cpp TBWTIndex.h BitRank.h Tools.h HuffWT.h \
  BlockArray.h
resultconvert.o: resultconvert.cpp
tconvert.o: tconvert.cpp
test.o: test.cpp NaiveForest.h Tools.h TBWTIndex.h BitRank.h HuffWT.h \
  BlockArray.h
tgenerate.o: tgenerate.cpp
traverse.o: traverse.cpp TBWTIndex.h BitRank.h Tools.h HuffWT.h \
  BlockArray.h
