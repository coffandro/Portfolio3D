#include <iostream>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <emscripten.h>

// Do something, if it goes wrong break the program
#define ASSERT(_e, ...) if (!(_e)) { fprintf(stderr, __VA_ARGS__); exit(1); }

// Define more readable definitions
typedef float    f32;
typedef double   f64;
typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int8_t   i8;
typedef int16_t  i16;
typedef int32_t  i32;
typedef int64_t  i64;
typedef size_t   usize;
typedef ssize_t  isize;

// Define window dimensions
#define WINDOW_WIDTH 768
#define WINDOW_HEIGHT 432
// Stolen from JDH
// Define vector 2 with float and int
typedef struct v2_s { f32 x, y; } v2;
typedef struct v2i_s { i32 x, y; } v2i;

// Get the dot product of two vectors
#define dot(v0, v1)                  \
    ({ const v2 _v0 = (v0), _v1 = (v1); (_v0.x * _v1.x) + (_v0.y * _v1.y); })

// Get the length of a vector???
#define length(v) ({ const v2 _v = (v); sqrtf(dot(_v, _v)); })

// normalises two vectors
#define normalize(u) ({              \
        const v2 _u = (u);           \
        const f32 l = length(_u);    \
        (v2) { _u.x / l, _u.y / l }; \
    })

// returns the smallest and biggest of the supplied values
#define min(a, b) ({ __typeof__(a) _a = (a), _b = (b); _a < _b ? _a : _b; })
#define max(a, b) ({ __typeof__(a) _a = (a), _b = (b); _a > _b ? _a : _b; })

// returns a -1 if a is negative, 0 if = 0 and 1 if positive
#define sign(a) ({                                       \
        __typeof__(a) _a = (a);                          \
        (__typeof__(a))(_a < 0 ? -1 : (_a > 0 ? 1 : 0)); \
    })

#define MAP_SIZE 10
static u8 MAPDATA[MAP_SIZE * MAP_SIZE] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 3, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    1, 0, 2, 0, 4, 4, 0, 0, 0, 1,
    1, 0, 0, 0, 4, 0, 0, 0, 0, 1,
    1, 0, 3, 0, 0, 0, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
};

// The struct containing most of the games actually... 
// important stuffs
struct {
    // Bare essentials to draw to SDL2
    SDL_Window *window;
    SDL_Texture *texture;
    SDL_Renderer *renderer;

    // a buffer to contain the pixels to draw on screen
    u32 pixels[WINDOW_WIDTH * WINDOW_HEIGHT];
    // a bool to keep track of whether to run the mainloop or not

    // three vectors to keep track of the position, direction
    // and something I am yet to understand.
    v2 pos, dir, plane;
} state;

// draw vertical line on x placement between y0 and y1
static void verline(int x, int y0, int y1, u32 color) {
    // for every int between y0 and y1 on the x cord place a
    // pixel on the screen with the color 
    for (int y = y0; y <= y1; y++) {
        state.pixels[(y * WINDOW_WIDTH) + x] = color;
    }
}

void setupSDL()
{
    // check if SDL can init video
    ASSERT(
        !SDL_Init(SDL_INIT_VIDEO),
        "SDL failed to initialize: %s\n",
        SDL_GetError());

    // create SDL window and check if working
    state.window =
        SDL_CreateWindow(
            "DEMO",
            SDL_WINDOWPOS_CENTERED_DISPLAY(0),
            SDL_WINDOWPOS_CENTERED_DISPLAY(0),
            WINDOW_WIDTH,
            WINDOW_HEIGHT,
            SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
    ASSERT(
        state.window,
        "failed to create SDL window: %s\n", SDL_GetError());

    // create SDL renderer and check if working
    state.renderer =
        SDL_CreateRenderer(state.window, -1, SDL_RENDERER_PRESENTVSYNC);
    ASSERT(
        state.renderer,
        "failed to create SDL renderer: %s\n", SDL_GetError());

    // create SDL texture and check if working
    state.texture =
        SDL_CreateTexture(
            state.renderer,
            SDL_PIXELFORMAT_ABGR8888,
            SDL_TEXTUREACCESS_STREAMING,
            WINDOW_WIDTH,
            WINDOW_HEIGHT);
    ASSERT(
        state.texture,
        "failed to create SDL texture: %s\n", SDL_GetError());
}

void init()
{
    // set start position of player, may consider removing
    // for 0,0 start
    state.pos = (v2) { 2, 2 };
    // set the current direction of the player
    state.dir = normalize(((v2) { -1.0f, 0.1f }));
    // still unsure
    state.plane = (v2) { 0.0f, 0.66f };
}

static void rotate(f32 rot) {
    const v2 d = state.dir, p = state.plane;
    state.dir.x = d.x * cos(rot) - d.y * sin(rot);
    state.dir.y = d.x * sin(rot) + d.y * cos(rot);
    state.plane.x = p.x * cos(rot) - p.y * sin(rot);
    state.plane.y = p.x * sin(rot) + p.y * cos(rot);
}

void process_input() {
    const f32
        rotspeed = 3.0f * 0.016f,
        movespeed = 3.0f * 0.016f;

    const u8 *keystate = SDL_GetKeyboardState(NULL);

    if (keystate[SDL_SCANCODE_Q] || keystate[SDL_SCANCODE_PAGEUP]) {
        rotate(+rotspeed);
    }

    if (keystate[SDL_SCANCODE_E] || keystate[SDL_SCANCODE_PAGEDOWN]) {
        rotate(-rotspeed);
    }

    if (keystate[SDL_SCANCODE_W] || keystate[SDL_SCANCODE_UP]) {
        state.pos.x += state.dir.x * movespeed;
        state.pos.y += state.dir.y * movespeed;
    }

    if (keystate[SDL_SCANCODE_S] || keystate[SDL_SCANCODE_DOWN]) {
        state.pos.x -= state.dir.x * movespeed;
        state.pos.y -= state.dir.y * movespeed;
    }

    if (keystate[SDL_SCANCODE_A] || keystate[SDL_SCANCODE_LEFT]) {
        state.pos.x += (state.dir.x * cos(1.5f) - state.dir.y * sin(1.5f)) * movespeed;
        state.pos.y += (state.dir.x * sin(1.5f) + state.dir.y * cos(1.5f)) * movespeed;
    }

    if (keystate[SDL_SCANCODE_D] || keystate[SDL_SCANCODE_RIGHT]) {
        state.pos.x += (state.dir.x * cos(-1.5f) - state.dir.y * sin(-1.5f)) * movespeed;
        state.pos.y += (state.dir.x * sin(-1.5f) + state.dir.y * cos(-1.5f)) * movespeed;
    }
}

void draw()
{
    // for every pixel on the x axis of the screen
    for (int x = 0; x < WINDOW_WIDTH; x++) {
        // x coordinate in space from [-1, 1]
        const f32 xcam = (2 * (x / (f32) (WINDOW_WIDTH))) - 1;

        // ray direction through this column
        const v2 dir = {
            state.dir.x + state.plane.x * xcam,
            state.dir.y + state.plane.y * xcam
        };

        // set local position and int position to state position
        v2 pos = state.pos;
        v2i ipos = { (int) pos.x, (int) pos.y };

        // distance ray must travel from one x/y side to the next
        const v2 deltadist = {
            static_cast<f32>(fabsf(dir.x) < 1e-20 ? 1e30 : fabsf(1.0f / dir.x)),
            static_cast<f32>(fabsf(dir.y) < 1e-20 ? 1e30 : fabsf(1.0f / dir.y)),
        };

        // distance from start position to first x/y side
        v2 sidedist = {
            deltadist.x * (dir.x < 0 ? (pos.x - ipos.x) : (ipos.x + 1 - pos.x)),
            deltadist.y * (dir.y < 0 ? (pos.y - ipos.y) : (ipos.y + 1 - pos.y)),
        };

        // integer step direction for x/y, calculated from overall diff
        const v2i step = { (int) sign(dir.x), (int) sign(dir.y) };

        // DDA hit
        struct { int val, side; v2 pos; } hit = { 0, 0, { 0.0f, 0.0f } };

        while (!hit.val) {
            if (sidedist.x < sidedist.y) {
                sidedist.x += deltadist.x;
                ipos.x += step.x;
                hit.side = 0;
            } else {
                sidedist.y += deltadist.y;
                ipos.y += step.y;
                hit.side = 1;
            }

            // kill player if they leave map
            //ASSERT(
            //    ipos.x >= 0
            //    && ipos.x < MAP_SIZE
            //    && ipos.y >= 0
            //    && ipos.y < MAP_SIZE,
            //    "DDA out of bounds");

            hit.val = MAPDATA[ipos.y * MAP_SIZE + ipos.x];
        }

        // set colors
        u32 color;
        switch (hit.val) {
        case 1: color = 0xFFFFFFFF; break;
        case 2: color = 0xFFFFFFFF; break;
        case 3: color = 0xFFFFFFFF; break;
        case 4: color = 0xFFFFFFFF; break;
        }

        // darken colors on y-sides
        if (hit.side == 1) {
            const u32
                br = ((color & 0xFF00FF) * 0xC0) >> 8,
                g  = ((color & 0x00FF00) * 0xC0) >> 8;

            color = 0xFF000000 | (br & 0xFF00FF) | (g & 0x00FF00);
        }

        // distance to hit
        const f32 dperp =
            hit.side == 0 ?
                (sidedist.x - deltadist.x)
                : (sidedist.y - deltadist.y);

        // perform perspective division, calculate line height relative to
        // screen center
        std::cout << std::to_string(dperp);

        const int
            h = (int) (WINDOW_HEIGHT / dperp),
            y0 = max((WINDOW_HEIGHT / 2) - (h / 2), 0),
            y1 = min((WINDOW_HEIGHT / 2) + (h / 2), WINDOW_HEIGHT - 1);

        verline(x, 0, y0, 0xFF202020);
        verline(x, y0, y1, color);
        verline(x, y1, WINDOW_HEIGHT - 1, 0xFF505050);
    }
}

void main_loop()
{
    process_input();

    memset(state.pixels, 0, sizeof(state.pixels));

    SDL_RenderClear(state.renderer);

    draw();

    SDL_UpdateTexture(state.texture, NULL, state.pixels, WINDOW_WIDTH * 4);
        SDL_RenderCopyEx(
            state.renderer,
            state.texture,
            NULL,
            NULL,
            0.0,
            NULL,
            SDL_FLIP_VERTICAL);
        SDL_RenderPresent(state.renderer);
}

void destroy()
{
    SDL_DestroyRenderer(state.renderer);
    SDL_DestroyWindow(state.window);
    SDL_Quit();
}

int main()
{
    setupSDL();
    init();

    emscripten_set_main_loop(main_loop, 60, true);

    destroy();
    return EXIT_SUCCESS;
}