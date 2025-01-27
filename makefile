BIN = bin
HTMLFILE = src/index.html

PATH_SDL_LINUX = /usr/include/SDL2

$(BIN):
	mkdir -p $@

web-release: $(BIN)
	rm -f $(BIN)/index.html
	cp -u $(HTMLFILE) $(BIN)
	emcc src/main.cpp \
	-DWEBBUILD=1 \
	-s WASM=1 \
	-s USE_SDL=2 \
	-s USE_SDL_IMAGE=2 \
	-s SDL2_IMAGE_FORMATS='["png"]' \
	-o bin/index.js

web-debug: $(BIN)
	rm -f $(BIN)/index.html
	emcc src/main.cpp \
	-DWEBBUILD=1 \
	-s WASM=1 \
	-s USE_SDL=2 \
	-s USE_SDL_IMAGE=2 \
	-s SDL2_IMAGE_FORMATS='["png"]' \
	-o bin/index.html

linux-release: $(BIN)
	gcc src/main.cpp \
	-lSDL2 \
	-lm \
	-lstdc++ \
	-o bin/main

all: $(BIN)	web-release	linux-release
