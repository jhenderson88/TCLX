CFLAGS = -std=c++0x -O2  -I$(PWD) -g `root-config --cflags`
SHAREDSWITCH = -shared -Wl,-soname,#NO ENDING SPACE
CPP = g++
COMPILESHARED = $(CPP) $(LFLAGS) $(SHAREDSWITCH)#NO ENDING SPACE

CFLAGS += -fPIC

OK_STRING="[OK]"
ERROR_STRING="[ERROR]"
WARN_STRING="[WARNING]"
COMP_STRING="Now Compiling "
FIN_STRING="Finished Building "

COM_COLOR=\033[0;34m
OBJ_COLOR=\033[0;36m
DICT_COLOR=\033[0;36m
OK_COLOR=\033[0;32m
ERROR_COLOR=\033[0;31m
WARN_COLOR=\033[0;33m
NO_COLOR=\033[m
FIN_COLOR=\033[3;34m
FIN_OBJ_COLOR=\033[3;32m

OBJECTS = TCLX.o TCLX_DataInput.o TReactions.o TExperiment.o TCLX_Dict.o 
INCDIR = $(PWD)
DICTOBJ = $(INCDIR)/TCLX.h $(INCDIR)/TCLX_DataInput.h $(INCDIR)/TReactions.h $(INCDIR)/TExperiment.h
DICTNAM = TCLX

all: libTExperiment.so
	@printf " ${FIN_COLOR}%s${FIN_OBJ_COLOR}%-30s ${NO_COLOR}\n" $(FIN_STRING) $^ ;

%.so: $(OBJECTS)
	@printf "${COM_COLOR}%s${OBJ_COLOR}%s${NO_COLOR}\n" $(COMP_STRING) "$@"
	@$(COMPILESHARED)$@ $(CFLAGS) -o$@ $(OBJECTS)

%.o: %.cxx
	@printf " ${COM_COLOR}%s ${OBJ_COLOR}%s${NO_COLOR}\n" $(COMP_STRING) $@ 
	@$(CXX) -c $^ $(CFLAGS) $(CPPFLAGS)

$(DICTNAM)_Dict.cxx: $(DICTOBJ)
	@printf "${COM_COLOR}%s${DICT_COLOR}%s${NO_COLOR}\n" $(COMP_STRING) $@
	@rootcint -f $@ -c -I$(INCDIR) -p $(DICTOBJ) linkdef.h

end: TExperiment
	@printf " ${WARN_COLOR}Compilation Success. woohoo!${NO_COLOR}\n\n"

clean:	
	$(RM) $(OBJECTS) *~ *Dict* *so

