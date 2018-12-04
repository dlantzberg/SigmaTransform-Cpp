#setup
VERSION = 3.02
CC      = g++
OBJ_DIR = obj
BIN_DIR = bin
OBJ1D 	= $(OBJ_DIR)/SigmaTransform1D.o $(OBJ_DIR)/SigmaTransform_util.o
OBJ2D 	= $(OBJ_DIR)/SigmaTransform2D.o $(OBJ_DIR)/SigmaTransform_util.o
ifdef OS
	#windows
	CFLAGS  = -std=gnu++11 -O3 -s -I"SigmaTransform/" -I"FFTW/"
	LDFLAGS = -lfftw3-3 -L"FFTW/"
	SYSTEM  = Windows
else
	#linux
	CFLAGS  = -std=gnu++11 -O3 -s -I"SigmaTransform/"
	LDFLAGS = -lfftw3 -lfftw3_threads -pthread
	SYSTEM  = Linux
endif
# targets
all: printSystem all1D all2D
	@echo "--- all done ---"
all1D: Example1D_STFT Example1D_ConstantQ Example1D_Wavelet Example1D_async Example1D_inline Example1D_threads
	@echo "--- done  1D ---"
all2D: Example2D_STFT Example2D_SIM2 Example2D_Curvelet Example2D_NPShearlet Example2D_Wavelet
	@echo "--- done  2D ---"
printSystem:
	@echo "--- OS: $(SYSTEM) ---"
# cleanup
.PHONY: clean
clean:
	rm $(OBJ_DIR)/*
	rm $(BIN_DIR)/*.bin
# link examples
Example1%: $(OBJ1D) $(OBJ_DIR)/Example1%.o
	$(CC) $(CFLAGS) -o $(BIN_DIR)/$@ $^ $(LDFLAGS)
Example2%: $(OBJ2D) $(OBJ_DIR)/Example2%.o
	$(CC) $(CFLAGS) -o $(BIN_DIR)/$@ $^ $(LDFLAGS)
# compile SigmaTransform-Code
$(OBJ_DIR)/%.o: SigmaTransform/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@
# compile examples
$(OBJ_DIR)/%.o: Examples/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@
