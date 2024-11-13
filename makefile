BIN = bin
HTMLFILE = src/index.html

$(BIN):
	mkdir -p $@

release: $(BIN)
	rm -f $(BIN)/index.html
	cp -u $(HTMLFILE) $(BIN)
	emcc src/main.cpp \
	-s WASM=1 \
	-s USE_SDL=2 \
	-s USE_SDL_IMAGE=2 \
	-s SDL2_IMAGE_FORMATS='["png"]' \
	-o bin/index.js

debug: $(BIN)
	rm -f $(BIN)/index.html
	emcc src/main.cpp \
	-s WASM=1 \
	-s USE_SDL=2 \
	-s USE_SDL_IMAGE=2 \
	-s SDL2_IMAGE_FORMATS='["png"]' \
	-o bin/index.html