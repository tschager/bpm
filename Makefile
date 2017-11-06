CFLAGS=-Wall -std=c++11
CFLAGSOPT=-O3
WITHMOD=-DMOD_TOLERANT
WITHMUT=-DMUT_TOLERANT

all: test main createIndex

main:
	$(CXX) $(CFLAGS) $(CFLAGSOPT) BPM.cpp src/*cpp -o BPM
	$(CXX) $(CFLAGS) $(CFLAGSOPT) $(WITHMOD) BPM.cpp src/*cpp -o BPM_Mod
	$(CXX) $(CFLAGS) $(CFLAGSOPT) $(WITHMUT) BPM.cpp src/*cpp -o BPM_Mut

test:
	$(CXX) $(CFLAGS) $(WITHMOD) $(WITHMUT) tests/unittests.cpp src/*cpp -o UnitTest

createIndex:
	$(CXX) $(CFLAGS) $(CFLAGSOPT) createDBIndex.cpp src/*cpp -o CreateIndex
	$(CXX) $(CFLAGS) $(CFLAGSOPT) $(WITHMOD) createDBIndex.cpp src/*cpp -o CreateIndex_Mod
	$(CXX) $(CFLAGS) $(CFLAGSOPT) $(WITHMOD) $(WITHMUT) createDBIndex.cpp src/*cpp -o CreateIndex_Mut

clean:
	 rm UnitTest BPM BPM_Mod BPM_Mut CreateIndex CreateIndex_Mod CreateIndex_Mut
