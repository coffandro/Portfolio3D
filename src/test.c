#include <stdio.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>

#include "assert.h"
#include "bitmap.h"
#include "rand.h"

// uncomment to enable lighting
// #define LIGHTMAPS

// number of lightmap pixels per 1 world unit
#define LIGHTMAP_SCALE 4

#define PI 3.14159265359f
#define TAU (2.0f * PI)
#define PI_2 (PI / 2.0f)
#define PI_4 (PI / 4.0f)

#define DEG2RAD(_d) ((_d) * (PI / 180.0f))
#define RAD2DEG(_d) ((_d) * (180f / PI))

#define SCREEN_WIDTH 384
#define SCREEN_HEIGHT 216

#define EYE_Z 2
#define HFOV DEG2RAD(90.0f)
#define VFOV 0.4f

#define ZNEAR 0.0001f
#define ZFAR  20.0f

// texture pixels per floor/ceiling unit
#define FLOOR_CEILING_PX_PER_UNIT 16

// texture pixels per sprite size unit
#define SPRITE_PX_PER_UNIT (1.0f / 16.0f)

#define UV_EPSILON 0.0001

#define LIGHT_MAX 255

typedef struct v2_s { f32 x, y; } v2;
typedef struct v3_s { f32 x, y, z; } v3;

typedef struct v2i_s { i32 x, y; } v2i;
typedef struct v3i_s { i32 x, y, z; } v3i;

struct Ray2 {
    v2 p, d;
};

struct Ray3 {
    v3 p, d;
};

// convert any v2* -> v2i
#define AS_V2I(_v) ({ TYPEOF(_v) __v = (_v); (v2i) { __v.x, __v.y }; })

// convert any v2* -> v2
#define AS_V2(_v) ({ TYPEOF(_v) __v = (_v); (v2) { __v.x, __v.y }; })

#define v2_eq(_a, _b, _eps) ({                                              \
        TYPEOF(_a) __a = (_a), __b = (_b);                                  \
        TYPEOF(_eps) __eps = (_eps);                                        \
        fabsf(__a.x - __b.x) < __eps && fabsf(__a.y - __b.y) < __eps;       \
    })

// 2D vector dot product
#define dot(_v0, _v1) ({                                                    \
    TYPEOF(_v0) __v0 = (_v0), __v1 = (_v1);                                 \
    (__v0.x * __v1.x) + (__v0.y * __v1.y); })

// 3D vector dot product
#define dot3(_v0, _v1) ({                                                   \
    TYPEOF(_v0) __v0 = (_v0), __v1 = (_v1);                                 \
    (__v0.x * __v1.x) + (__v0.y * __v1.y) + (__v0.z * __v1.z); })

// 2D vector length
#define length(_vl) ({ TYPEOF(_vl) __vl = (_vl); sqrtf(dot(__vl, __vl)); })

// 3D vector length
#define length3(_vl) ({ TYPEOF(_vl) __vl = (_vl); sqrtf(dot3(__vl, __vl)); })

// 2D vector normalization
#define normalize(_vn) ({                                                   \
        TYPEOF(_vn) __vn = (_vn);                                           \
        const f32 l = length(__vn);                                         \
        (TYPEOF(_vn)) { __vn.x / l, __vn.y / l }; })

// 3D vector normalization
#define normalize3(_vn) ({                                                  \
        TYPEOF(_vn) __vn = (_vn);                                           \
        const f32 l = length3(__vn);                                        \
        (TYPEOF(_vn)) { __vn.x / l, __vn.y / l, __vn.z / l }; })

// 2D vector cross product
#define cross(_vc0, _vc1) ({                                                \
        TYPEOF(_vc0) __vc0 = (_vc0), __vc1 = (_vc1);                        \
        (__vc0.x * __vc1.y) - (__vc1.x * __vc0.y); })

// 2D vector rotation
#define rotate(_vr, _th) ({                                                 \
        TYPEOF(_vr) __vr = (_vr);                                           \
        TYPEOF(_th) __th = (_th);                                           \
        (v2) {                                                              \
            (__vr.x * cos(__th)) - (__vr.y * fastsin(__th)),                \
            (__vr.x * fastsin(__th)) + (__vr.y * cos(__th)),                \
        };                                                                  \
    })

// see: https://en.wikipedia.org/wiki/Line–line_intersection

// intersect two infinite lines
#define intersect_lines(_a0, _a1, _b0, _b1) ({                              \
        TYPEOF(_a0) _l00 = (_a0), _l01 = (_a1), _l10 = (_b0), _l11 = (_b1); \
        TYPEOF(_l00.x) _d =                                                 \
            ((_l00.x - _l01.x) * (_l10.y - _l11.y))                         \
                - ((_l00.y - _l01.y) * (_l10.x - _l11.x)),                  \
            _xl0 = cross(_l00, _l01),                                       \
            _xl1 = cross(_l10, _l11);                                       \
        (TYPEOF(_a0)) {                                                     \
            ((_xl0 * (_l10.x - _l11.x)) - ((_l00.x - _l01.x) * _xl1)) / _d, \
            ((_xl0 * (_l10.y - _l11.y)) - ((_l00.y - _l01.y) * _xl1)) / _d  \
        };                                                                  \
    })

// TODO: optimize to fail early
// intersect two line segments, returns (nan, nan) if no intersection exists
#define intersect_segs(_a0, _a1, _b0, _b1) ({                               \
        TYPEOF(_a0) _l00 = (_a0), _l01 = (_a1), _l10 = (_b0), _l11 = (_b1); \
        TYPEOF(_l00.x)                                                      \
            _d = ((_l00.x - _l01.x) * (_l10.y - _l11.y))                    \
                    - ((_l00.y - _l01.y) * (_l10.x - _l11.x)),              \
            _t =                                                            \
                (((_l00.x - _l10.x) * (_l10.y - _l11.y))                    \
                    - ((_l00.y - _l10.y) * (_l10.x - _l11.x)))              \
                    / _d,                                                   \
            _u =                                                            \
                (((_l00.x - _l10.x) * (_l00.y - _l01.y))                    \
                    - ((_l00.y - _l10.y) * (_l00.x - _l01.x)))              \
                    / _d;                                                   \
        (fabsf(_d) < 0.00001f) ?                                            \
            ((v2) { NAN, NAN })                                             \
            : ((_t >= 0 && _t <= 1 && _u >= 0 && _u <= 1) ?                 \
                ((v2) {                                                     \
                    _l00.x + (_t * (_l01.x - _l00.x)),                      \
                    _l00.y + (_t * (_l01.y - _l00.y)) })                    \
                : ((v2) { NAN, NAN }));                                     \
    })

// TODO: optimize to fail early
// true if _p is in triangle _a, _b, _c
#define point_triangle(_p, _a, _b, _c) ({                                   \
        TYPEOF(_p) __p = (_p), __a = (_a), __b = (_b), __c = (_c);          \
        const f32                                                           \
            _d = ((__b.y - __c.y) * (__a.x - __c.x)                         \
                  + (__c.x - __b.x) * (__a.y - __c.y)),                     \
            _x = ((__b.y - __c.y) * (__p.x - __c.x)                         \
                  + (__c.x - __b.x) * (__p.y - __c.y)) / _d,                \
            _y = ((__c.y - __a.y) * (__p.x - __c.x)                         \
                  + (__a.x - __c.x) * (__p.y - __c.y)) / _d,                \
            _z = 1 - _x - _y;                                               \
        (_x > 0) && (_y > 0) && (_z > 0);                                   \
    })

// -1 right, 0 on, 1 left
#define point_side(_p, _a, _b) ({                                           \
        TYPEOF(_p) __p = (_p), __a = (_a), __b = (_b);                      \
        -(((__p.x - __a.x) * (__b.y - __a.y))                               \
            - ((__p.y - __a.y) * (__b.x - __a.x)));                         \
    })

// returns fractional part of float _f
#define fract(_f) ({ TYPEOF(_f) __f = (_f); __f - ((i64) (__f)); })

// -1, 0, 1 depending on _f < 0, _f == 0, _f > 0
#define sign(_f) ({                                                         \
    TYPEOF(_f) __f = (_f);                                                  \
    (TYPEOF(_f)) (__f < 0 ? -1 : (__f > 0 ? 1 : 0)); })

#define min(_a, _b) ({                                                      \
        TYPEOF(_a) __a = (_a), __b = (_b);                                  \
        __a < __b ? __a : __b; })

#define max(_a, _b) ({                                                      \
        TYPEOF(_a) __a = (_a), __b = (_b);                                  \
        __a > __b ? __a : __b; })

// clamp _x such that _x is in [_mi.._ma]
#define clamp(_x, _mi, _ma) (min(max(_x, _mi), _ma))

// swap _a and _b
#define swap(_a, _b) ({ TYPEOF(_a) _x = (_a); _a = (_b); _b = (_x); })

// if _x is nan, returns _alt otherwise returns _x
#define ifnan(_x, _alt) ({ TYPEOF(_x) __x = (_x); isnan(__x) ? (_alt) : __x; })

// lerp from _a -> _b by _t
#define lerp(_a, _b, _t) ({                                       \
        TYPEOF(_t) __t = (_t);                                    \
        (TYPEOF(_a)) (((_a) * (1 - __t)) + ((_b) * __t));         \
    })

// multiple B, G, R of _abgr with _x and clamp in [0x00..0xFF]
#define abgr_mul(_abgr, _x) ({                                    \
        TYPEOF(_abgr) __abgr = (_abgr);                           \
        (__abgr & 0xFF000000)                                     \
        | ((((int) (((__abgr >> 16) & 0xFF) * _x)) & 0xFF) << 16) \
        | ((((int) (((__abgr >> 8) & 0xFF) * _x)) & 0xFF) << 8)   \
        | ((((int) (((__abgr >> 0) & 0xFF) * _x)) & 0xFF) << 0);  \
    })

// blend two abgr colors with _c0 as "dst" and _c1 as "src"
#define abgr_blend(_c0, _c1) ({                                     \
        TYPEOF(_c0) __c0 = (_c0), __c1 = (_c1);                     \
        const f32 _t = (__c1 >> 24) / 255.0f;                       \
        (__c1 >> 24)                                                \
        | (lerp((__c0 >> 16) & 0xFF, (__c1 >> 16) & 0xFF) << 16, _t)\
        | (lerp((__c0 >> 8) & 0xFF, (__c1 >> 8) & 0xFF) << 8, _t)   \
        | (lerp((__c0 >> 0) & 0xFF, (__c1 >> 0) & 0xFF) << 0, _t)   \
    })

ALWAYS_INLINE void memset32(void *dst, u32 val, usize n) {
    u16 *dst32 = dst;
    while (n--) { *dst32++ = val; }
}

ALWAYS_INLINE void memset16(void *dst, u16 val, usize n) {
    u16 *dst16 = dst;
    while (n--) { *dst16++ = val; }
}

ALWAYS_INLINE void memsetf64(void *dst, f64 val, usize n) {
    f64 *dst64 = dst;
    while (n--) { *dst64++ = val; }
}

ALWAYS_INLINE void memsetf32(void *dst, f32 val, usize n) {
    f32 *dst32 = dst;
    while (n--) { *dst32++ = val; }
}

struct Light {
    int sector;
    v3 pos;
    f32 power;
};

struct Wall {
    // READ VALUES
    v2 a, b;            // points of wall, front face is left of a -> b
    int portal;         // if not SECTOR_NONE, sector wall is portal to
    int material;       // index of material

    // COMPUTED VALUES
    int index;          // index in walls array
    int sector;         // sector this wall belongs to
    v2 normal;          // normal (left of a -> b)
    v2 d;               // b - a
    f32 inv_d_dot_d;    // 1.0 / (d dot d), used for fast point-nearest-to-wall
    f32 len;            // length(d)
    f32 height;         // sector ceil - sector floor TODO REMOVE
    int portal_wall;    // wall to which this wall is a portal
    bool portal_same;   // if true, the wall this wall is a portal do does NOT
                        // have the same location as this wall

#ifdef LIGHTMAPS
    struct {
        u8 *data;
        v2i size;
    } lightmap;
#endif
};

enum PlaneType {
    PLANE_TYPE_FLOOR,
    PLANE_TYPE_CEIL
};

// floor or ceiling plane
struct Plane {
    // READ VALUES
    f32 z;                  // plane offset from 0
    int texture;            // texture
    v2 scale, offset;       // texture scale and offset

    // COMPUTED VALUES
    enum PlaneType type;    // type of plane

#ifdef LIGHTMAPS
    struct {
        u8 *data;
        v2i size;
    } lightmap;
#endif
};

// sector id for "no sector"
#define SECTOR_NONE 0
#define SECTOR_FIRST 1

// maximum number of sectors in a level
#define SECTOR_MAX 128

// maximum number of allowed portals to one sector
#define PORTAL_MAX 16

struct Sector {
    // READ VALUES
    int id;                    // sector id (>0 if valid)
    usize first_wall, n_walls; // first wall index, number of walls
    u8 light;                  // sector light value
    struct Plane floor, ceil;  // floor/ceiling planes

    // COMPUTED VALUES
    v2 min, max, span; // minimum/maximum bounds of this sector and max - min
};

enum TextureId {
    TEXTURE_NONE = 0,
    TEXTURE_TEST0,
    TEXTURE_TEST1,
    TEXTURE_TEST2,
    TEXTURE_TEST3,
    TEXTURE_TEST4,
    TEXTURE_COUNT
};

struct Texture {
    int id;
    u32 *data;
    v2i offset, size;
    usize pitch;
};

#define MATERIAL_NONE 0
struct Material {
    int id, texture_id;
    v2 scale;
};

#define TRIG_TAB_SIZE 2048
static f32 sintab[TRIG_TAB_SIZE], tantab[TRIG_TAB_SIZE];

// portal clipping data from one sector (from) to another
struct PortalClip {
    int from;
    int x0, x1, yf0, yc0;
    f32 mf, mc;
};

// clipping data for a whole sector with up to PORTAL_MAX portals
struct SectorClip {
    struct PortalClip arr[PORTAL_MAX];  // clipping data for each portal
    usize n;                            // length of arr
    int x0, x1;                         // x0 -> x1 of entire clipping area
};

static struct {
    SDL_Window *window;
    SDL_Renderer *renderer;
    SDL_Texture *texture, *debug;
    u32 *pixels, *texture_pixels;
    bool quit;

    struct { struct Sector arr[32]; usize n; } sectors;
    struct { struct Wall arr[128]; usize n; } walls;
    struct { struct Material arr[128]; usize n; } materials;
    struct Texture textures[TEXTURE_COUNT];

    f32 wall_dist[SCREEN_WIDTH];
    u16 y_lo[SCREEN_WIDTH], y_hi[SCREEN_WIDTH];

    // bitmap where entry for [SECTOR ID] is true if that sector was drawn at
    // all this frame
    u8 sectdraw[BITMAP_SIZE_TO_BYTES(SECTOR_MAX)];

    // clipping data for each (drawn) sector
    struct SectorClip sectclip[SECTOR_MAX];

    // 2D bitmap indexed [SECTOR ID A][SECTOR ID B], each entry indicating that
    // a portal from A -> B has been drawn. used to prevent the case that a
    // portal from A -> B has been drawn and B is trying to draw a portal back
    // to A (causes infinite render times :( )
    u8 sectport[BITMAP_SIZE_TO_BYTES(SECTOR_MAX * SECTOR_MAX)];

    struct {
        v2 pos;
        f32 z;
        f32 angle, anglecos, anglesin;
        int sector;
    } camera;

    bool sleepy;
} state;

static struct Light g_lights[1] = {
    {
        5,
        { 4.06, 4.04, 2.04 },
        12.0
    }
};

static void make_trigtabs() {
    for (usize i = 0; i < TRIG_TAB_SIZE; i++) {
        // 0..TAU
        const f32 a = (i * TAU) / (f32) (TRIG_TAB_SIZE);
        sintab[i] = sin(a);
        tantab[i] = tan(a);
    }
}

ALWAYS_INLINE f32 fastsin(f32 a) {
    const long i = lround(((a) * (TRIG_TAB_SIZE / 2)) / PI);
    if (i < 0) {
        return sintab[(TRIG_TAB_SIZE - ((-i) & (TRIG_TAB_SIZE - 1))) & (TRIG_TAB_SIZE - 1)];
    } else {
        return sintab[i & (TRIG_TAB_SIZE - 1)];
    }
}


ALWAYS_INLINE f32 fastcos(f32 a) {
    const long i = lround(((a) * (TRIG_TAB_SIZE / 2)) / PI);
    if (i < 0) {
        return sintab[((-i) + (TRIG_TAB_SIZE / 4)) & (TRIG_TAB_SIZE - 1)];
    } else {
        return sintab[(i + (TRIG_TAB_SIZE / 4)) & (TRIG_TAB_SIZE - 1)];
    }
}

ALWAYS_INLINE f32 fasttan(f32 a) {
    const long i = lround(((a) * (TRIG_TAB_SIZE / 2)) / PI);
    if (i < 0) {
        return tantab[(TRIG_TAB_SIZE - ((-i) & (TRIG_TAB_SIZE - 1))) & (TRIG_TAB_SIZE - 1)];
    } else {
        return tantab[i & (TRIG_TAB_SIZE - 1)];
    }
}

// optimal 3-term polynomial approximation of atan(), error is max 0.0008
// stackoverflow.com/questions/42537957
ALWAYS_INLINE f32 fastatan(f32 x) {
#define FT_A 0.0776509570923569
#define FT_B -0.287434475393028
#define FT_C (PI_4 - FT_A - FT_B)
    const f32 xx = x * x;
    return ((FT_A * xx + FT_B) * xx + FT_C) * x;
#undef FT_A
#undef FT_B
#undef FT_C
}

// find point on wall nearest to p
// see: math.stackexchange.com/questions/2193720
ALWAYS_INLINE v2 wall_point_nearest_to(const struct Wall *wall, v2 p) {
    const f32 t =
        -dot(wall->d, ((v2) { wall->a.x - p.x, wall->a.y - p.y }))
            * wall->inv_d_dot_d;

    if (t >= 0 && t <= 1) {
        return (v2) { wall->a.x + t * wall->d.x, wall->a.y + t * wall->d.y };
    }

    const f32
        l_a = length(((v2) { wall->a.x - p.x, wall->a.y - p.y })),
        l_b = length(((v2) { wall->b.x - p.x, wall->b.y - p.y }));

    return l_a < l_b ? AS_V2(wall->a) : AS_V2(wall->b);
}

// convert angle in [-(HFOV / 2)..+(HFOV / 2)] to X coordinate relative to
// screen center
ALWAYS_INLINE int screen_angle_to_unbounded_x(f32 angle) {
    // convert to [-PI/4..+PI/4]
    angle = (((angle + (HFOV / 2.0f)) / HFOV) * PI_2) - PI_4;

    // use -tan(angle) since screen x coordinates go the opposite direction of
    // relative camera angle (i.e. SCREEN_WIDTH - 1 is at -HFOV/2)
    // (NOLINTNEXTLINE)
    return (int) ((SCREEN_WIDTH / 2) * (1.0f + -fasttan(angle)));
}

ALWAYS_INLINE int screen_angle_to_x(f32 angle) {
    return clamp(screen_angle_to_unbounded_x(angle), 0, SCREEN_WIDTH - 1);
}

// convert screen X to [-(HFOV / 2)..+(HFOV / 2)] (inverse of angle -> x)
ALWAYS_INLINE f32 x_to_screen_angle(int x) {
    // (NOLINTNEXTLINE) convert back to tan result
    const f32 a = -((x / (f32) (SCREEN_WIDTH / 2)) - 1.0f);

    // arctan [-PI/4..PI/4] -> [-1, 1] -> [-HFOV/2..+HFOV/2]
    return (((fastatan(a) + PI_4) / PI_2) * HFOV) - (HFOV / 2.0f);
}

// noramlize angle to +/-PI
ALWAYS_INLINE f32 normalize_angle(f32 a) {
    return a - (TAU * floorf((a + PI) / TAU));
}

// camera space -> world space (un-rotate and un-translate)
ALWAYS_INLINE v2 camera_pos_to_world(v2 p) {
    const v2
        pr = {
            p.x * state.camera.anglesin + p.y * state.camera.anglecos,
            p.y * state.camera.anglesin - p.x * state.camera.anglecos,
        };
    return (v2) { pr.x + state.camera.pos.x, pr.y + state.camera.pos.y };
}

// world space -> camera space (translate and rotate)
// rotates points by the matrix
// | sin(a) -cos(a) |
// | cos(a)  sin(a) |
// such that y = 0 becomes x = 0 (as if the world was relative to the y-axis)
ALWAYS_INLINE v2 world_pos_to_camera(v2 p) {
    const v2 u = { p.x - state.camera.pos.x, p.y - state.camera.pos.y };
    return (v2) {
        u.x * state.camera.anglesin - u.y * state.camera.anglecos,
        u.y * state.camera.anglesin + u.x * state.camera.anglecos,
    };
}

// convert a wall position AND a screen "y" (floor or ceiling) to a world point
// based on the floor or ceiling height "z"
//
// HOW THIS WORKS
//
// to find world point p = (p_x, p_y)
//
// let
//   v ("dir_to_wall") be the direction vector from the camera to the wall slice
//   w ("wall_point") be the camera-space point of the wall slice
//   y be some screen y-coordinate corresponding to the wall point
//   s_y (sy0, sy1 in main rendering function) by the y-scale for this wall
//   slice
//   z be the ceiling or floor height of the y input
//
// then:
//
// y = (HEIGHT / 2) + (z - Z_EYE) * s_y (from rendering calculations)
//
// s_y = (y - (HEIGHT / 2)) / (z - Z_EYE)
// s_y = (VFOV * HEIGHT) / p_y
//
// p_y = (VFOV * HEIGHT) / s_y
// p_y = ((VFOV * HEIGHT) * (z - Z_EYE)) / (y - (HEIGHT / 2))
//
// now, we can take advantage of the fact that the vertical slices of floor/
// ceiling which are drawn travel ALONG the vector -v (wall to camera) and say:
//
// w - t * v = p
// for some scaling factor "t"
//
// expanding this out component wise, we have
// w_x - t * v_x = p_x
// w_y - t * v_y = p_y
//
// since we have solved for p_y above, we can use this to solve for p_x
//
// p_x = w_x - ((w_y - p_y) / v_y) * v_x
//
// giving us p = (p_x, p_y) in camera space
//
// NOTE: this is (I think) not the fastest way to do this as it adds an extra
// two (!) divisions per-pixel :(
ALWAYS_INLINE v2 wall_floor_ceiling_to_camera_space(
    v2 dir_to_wall,
    v2 wall_point,
    int y,
    f32 z) {
    const int
        d = (y - (SCREEN_HEIGHT / 2)),
        e = d == 0 ? 1 : d;

    const f32 py = ((VFOV * SCREEN_HEIGHT) * (z - state.camera.z)) / ((f32) e);
    return (v2) {
        wall_point.x - (((wall_point.y - py) / dir_to_wall.y) * dir_to_wall.x),
        py
    };
}

static void present();

// load sectors from file -> state
static int load_sectors(const char *path) {
    // sector 0, material 0 do not exist
    state.sectors.n = 1;
    state.materials.n = 1;

    FILE *f = fopen(path, "r");
    if (!f) { return -1; }

    int retval = 0;
    enum { SCAN_SECTOR, SCAN_WALL, SCAN_NONE } ss = SCAN_NONE;

    char line[1024], buf[64];
    while (fgets(line, sizeof(line), f)) {
        const char *p = line;
        while (isspace(*p)) {
            p++;
        }

        // skip line, empty or comment
        if (!*p || *p == '#') {
            continue;
        } else if (*p == '[') {
            strncpy(buf, p + 1, sizeof(buf));
            const char *section = strtok(buf, "]");
            if (!section) { retval = -2; goto done; }

            if (!strcmp(section, "SECTOR")) { ss = SCAN_SECTOR; }
            else if (!strcmp(section, "WALL")) { ss = SCAN_WALL; }
            else { retval = -3; goto done; }
        } else {
            switch (ss) {
            case SCAN_WALL: {
                struct Wall *wall = &state.walls.arr[state.walls.n++];
                if (sscanf(
                        p,
                        "%f %f %f %f %d",
                        &wall->a.x,
                        &wall->a.y,
                        &wall->b.x,
                        &wall->b.y,
                        &wall->portal)
                        != 5) {
                    retval = -4; goto done;
                }
                wall->material = 1;
            }; break;
            case SCAN_SECTOR: {
                struct Sector *sector = &state.sectors.arr[state.sectors.n++];
                if (sscanf(
                        p,
                        "%d %" PRIusize " %" PRIusize " %f %f %" PRIu8,
                        &sector->id,
                        &sector->first_wall,
                        &sector->n_walls,
                        &sector->floor.z,
                        &sector->ceil.z,
                        &sector->light)
                        != 6) {
                    retval = -5; goto done;
                }

                sector->floor.texture = (TEXTURE_TEST0 + ((sector->id - 1) % 4));
                sector->floor.scale = (v2) { 1, 1 };
                sector->floor.offset = (v2) { 0, 0 };

                sector->ceil.texture = 3;
                sector->ceil.scale = (v2) { 1, 1 };
                sector->ceil.offset = (v2) { 0, 0 };
            }; break;
            default: retval = -6; goto done;
            }
        }
    }

    if (ferror(f)) { retval = -128; goto done; }

    // precalculate sector data
    for (usize i = SECTOR_FIRST; i < state.sectors.n; i++) {
        struct Sector *sector = &state.sectors.arr[i];

        sector->floor.type = PLANE_TYPE_FLOOR;
        sector->ceil.type = PLANE_TYPE_CEIL;

        v2 p_min = { 1e10, 1e10 }, p_max = { 0, 0 };

        for (usize j = 0; j < sector->n_walls; j++) {
            struct Wall *wall = &state.walls.arr[sector->first_wall + j];
            wall->index = sector->first_wall + j;
            wall->sector = sector->id;

            p_min.x = min(p_min.x, min(wall->a.x, wall->b.x));
            p_min.y = min(p_min.y, min(wall->a.y, wall->b.y));
            p_max.x = max(p_max.x, max(wall->a.x, wall->b.x));
            p_max.y = max(p_max.y, max(wall->a.y, wall->b.y));

            // compute world-space normal of front side (right of a -> b)
            // use the fact that for a line segment (x0, y0) -> (x1, y1), with
            // dx = x1 - x0, dy = y1 - y0, then its 2D normals are (-dy, dx) and
            // (dy, -dx). all walls are perpendicular so the "z" component of
            wall->normal =
                normalize(
                    ((v2) {
                         wall->b.y - wall->a.y,
                        -(wall->b.x - wall->a.x)
                    }));
            wall->d =
                (v2) {
                    wall->b.x - wall->a.x,
                    wall->b.y - wall->a.y
                };
            wall->inv_d_dot_d = 1.0f / dot(wall->d, wall->d);
            wall->len = length(wall->d);
            wall->height = sector->ceil.z - sector->floor.z;

            if (wall->portal) {
                // find corresponding wall in portal sector
                const struct Sector *neighbor =
                    &state.sectors.arr[wall->portal];

                const struct Wall *portal_wall = NULL;

                for (usize k = 0; k < neighbor->n_walls; k++) {
                    const struct Wall *w  =
                        &state.walls.arr[neighbor->first_wall + k];

                    if (w->portal == sector->id) {
                        portal_wall = w;
                        break;
                    }
                }

                if (!portal_wall) {
                    WARN(
                        "wall %d has portal to sector %d,"
                        " but sector has no portal back!",
                        sector->first_wall + j,
                        wall->portal);
                }

                wall->portal_wall = portal_wall->index;

                // set portal_same if the portal wall and this wall are in the
                // same position
                wall->portal_same =
                    v2_eq(wall->a, portal_wall->a, 0.00001f)
                    && v2_eq(wall->b, portal_wall->b, 0.00001f);
            }
        }

        sector->min = p_min;
        sector->max = p_max;
        sector->span = (v2) { p_max.x - p_min.x, p_max.y - p_min.y };
    }
done:
    fclose(f);
    return retval;
}

// clamp point to be inside of sector
static v2 sector_clamp_point(const struct Sector *sector, v2 p) {
    bool in_sector = true;

    usize n = 0;
    v2 qs[sector->n_walls];

    for (usize i = 0; i < sector->n_walls; i++) {
        const struct Wall *wall = &state.walls.arr[sector->first_wall + i];

        // if point is on the left side of line, add candidate clamped to wall
        if (point_side(p, AS_V2(wall->a), AS_V2(wall->b)) > 0) {
            in_sector = false;
            qs[n++] = wall_point_nearest_to(wall, p);
        }
    }

    // no need to clamp
    if (in_sector) {
        return p;
    }

    // find nearest point in qs to p
    v2 r = qs[0];
    f32 d_r = 1e30f;

    for (usize i = 0; i < n; i++) {
        const f32 d = length(((v2) { qs[i].x - p.x, qs[i].y - p.y }));
        if (d < d_r) {
            r = qs[i];
            d_r = d;
        }
    }

    return r;
}

// clamp 3D point to be inside of sector
static v3 sector_clamp_point3(const struct Sector *sector, v3 p) {
    const v2 q = sector_clamp_point(sector, (v2) { p.x, p.y });
    return (v3) { q.x, q.y, clamp(p.z, sector->floor.z, sector->ceil.z) };
}

// point is in sector if it is on the right side of all walls
static bool point_in_sector(const struct Sector *sector, v2 p) {
    for (usize i = 0; i < sector->n_walls; i++) {
        const struct Wall *wall = &state.walls.arr[sector->first_wall + i];

        if (point_side(p, AS_V2(wall->a), AS_V2(wall->b)) > 0) {
            return false;
        }
    }

    return true;
}

struct RaySectorIntersection {
    const struct Wall *wall;
    v2 p;
};

bool intersect_ray_sector(
    struct Ray2 ray,
    const struct Sector *sector,
    struct RaySectorIntersection *out) {
    const f32 l = max(sector->span.x * 4, sector->span.y * 4);
    const v2 p1 = { ray.p.x + (ray.d.x * l), ray.p.y + (ray.d.y * l) };

    for (usize i = 0; i < sector->n_walls; i++) {
        const struct Wall *wall = &state.walls.arr[sector->first_wall + i];

        const v2 t =
            intersect_segs(
                AS_V2(wall->a),
                AS_V2(wall->b),
                ray.p,
                p1);

        if (!isnan(t.x)) {
            *out = (struct RaySectorIntersection) { wall, t };
            return true;
        }
    }

    return false;
}

struct RayLevelIntersection {
    const struct Sector *sector;
    const struct Wall *wall;
    v3 p;
    f32 dist;
};

bool intersect_ray_level(
    struct Ray3 ray,
    int sector_ray,
    int stop_at_sector,
    struct RayLevelIntersection *out) {
#define RAYCAST_EPSILON     0.0001f
#define RAYCAST_MAX_SECTORS 16

    // current ray sector
    const struct Sector *sector = &state.sectors.arr[sector_ray];

    // don't attempt raycast if in stopping sector
    if (sector_ray == stop_at_sector) {
        *out = (struct RayLevelIntersection) { sector, NULL, ray.p, 0.0f };
        return true;
    }

    // current ray length
    f32 dist = 0.0f;

    // current ray point, start by sliding along direction a little bit so that
    // points at sector edges don't collide with themselves
    // also clamp to sector, ray may have left after epsilon adjustment
    v3 p =
        sector_clamp_point3(
            sector,
            (v3) {
                ray.p.x + (ray.d.x * RAYCAST_EPSILON),
                ray.p.y + (ray.d.y * RAYCAST_EPSILON),
                ray.p.z + (ray.d.z * RAYCAST_EPSILON)
            });

    // number of sector hits
    usize count = 0;

    struct RaySectorIntersection rsi;
    while (count < RAYCAST_MAX_SECTORS
           && intersect_ray_sector(
               (struct Ray2) { { p.x, p.y }, { ray.d.x, ray.d.y } },
               sector,
               &rsi)) {
        const struct Sector *neighbor =
            rsi.wall->portal ?
                &state.sectors.arr[rsi.wall->portal]
                : NULL;

        // find ray hit "t" (where hit = t * ray_direction)
        // choose whichever direction travels furthest (x or y) to avoid 0-div
        // issues
        f32 t =
            fabsf(ray.d.x) > fabsf(ray.d.y) ?
                ((rsi.p.x - p.x) / ray.d.x)
                : ((rsi.p.y - p.y) / ray.d.y);

        // find xy/z distance to hit, z of hit
        const f32 z_hit = p.z + (t * ray.d.z);
        f32 z = z_hit;

        if (z_hit < sector->floor.z) {
            // hit on floor
            t = (sector->floor.z - p.z) / ray.d.z;
            z = sector->floor.z;
        } else if (z_hit > sector->ceil.z) {
            // hit on ceiling
            t = (sector->ceil.z - p.z) / ray.d.z;
            z = sector->ceil.z;
        } else if (neighbor && z_hit < neighbor->floor.z) {
            // hit on lower portal wall
            t = (neighbor->floor.z - p.z) / ray.d.z;
            z = p.z + (t * ray.d.z);
        } else if (neighbor && z_hit > neighbor->ceil.z) {
            // hit on upper portal wall
            t = (neighbor->ceil.z - p.z) / ray.d.z;
            z = p.z + (t * ray.d.z);
        }

        p = (v3) { rsi.p.x, rsi.p.y, z };
        dist += length3(((v3) { t * ray.d.x, t * ray.d.y, t * ray.d.z }));

        if (z_hit < sector->floor.z || z_hit > sector->ceil.z) {
            // hit on floor or ceiling
            *out = (struct RayLevelIntersection) { sector, NULL, p, dist };
            return true;
        }

        if (!rsi.wall->portal
            || z_hit < neighbor->floor.z
            || z_hit > neighbor->ceil.z) {
            // hit solid wall or above/below portal wall
            *out = (struct RayLevelIntersection) { sector, rsi.wall, p, dist };
            return true;
        }

        // should we stop at this neighbor?
        if (neighbor->id == stop_at_sector) {
            *out = (struct RayLevelIntersection) { neighbor, NULL, p, dist };
            return true;
        }

        v2 xy = { p.x, p.y };

        // if wall position is not exactly the same we need to move the ray in
        // the xy-plane
        if (!rsi.wall->portal_same) {
            const struct Wall *neighbor_wall =
                &state.walls.arr[rsi.wall->portal_wall];

            // find our hit "u" on this wall such that for this wall A -> B
            // hit = A + u * (B - A)
            // we can use the same u on the other wall C -> D in the form
            // hit_other = C + u * (D - C)
            // in case the two walls are not in the same position despite being
            // portals
            const f32 u =
                length(
                    ((v2) {
                        rsi.p.x - rsi.wall->a.x,
                        rsi.p.y - rsi.wall->a.y }))
                    / rsi.wall->len;

            // and use it to get the point on the neighbor wall as they
            // *might* not be the same
            xy = (v2) {
                neighbor_wall->a.x + u * neighbor_wall->d.x,
                neighbor_wall->a.y + u * neighbor_wall->d.y
            };
        }

        // use same epsilon trick to avoid colliding with portal again
        p = (v3) {
            p.x + (xy.x * RAYCAST_EPSILON),
            p.y + (xy.y * RAYCAST_EPSILON),
            p.z + (ray.d.z * RAYCAST_EPSILON)
        };

        // continue raycast from p_neighbor
        sector = neighbor;

        count++;
    }

    *out = (struct RayLevelIntersection) { NULL, NULL, { NAN, NAN, NAN }, NAN };
    return false;

#undef RAYCAST_EPSILON
#undef RAYCAST_MAX_SECTORS
}

#ifdef LIGHTMAPS
static u8 calculate_point_light(
    const struct Light *lights,
    usize n_lights,
    int sector_point,
    v3 p,
    v3 n) {
    u8 l = 0;

    for (usize i = 0; i < n_lights; i++) {
        const struct Light *light = &lights[i];

        const v3 p_to_l = {
            light->pos.x - p.x,
            light->pos.y - p.y,
            light->pos.z - p.z
        };

        if (length(p_to_l) > light->power) {
            continue;
        }

        const v3 d_to_light = normalize3(p_to_l);

        struct RayLevelIntersection rli;
        const bool hit =
            intersect_ray_level(
                (struct Ray3) {
                    p,
                    d_to_light
                },
                sector_point,
                light->sector,
                &rli);

        if (!hit || (hit && rli.sector->id != light->sector)) {
            continue;
        }

        // add last distance from point to light
        const f32 dist =
            rli.dist
                + length3(((v3) {
                    light->pos.x - rli.p.x,
                    light->pos.y - rli.p.y,
                    light->pos.z - rli.p.z,
                }));

        const f32
            d = 1.0f - clamp(dist / light->power, 0.0f, 1.0f),
            n_dot_l = clamp(dot3(d_to_light, n), 0.0f, 1.0f);

        l = min(l + (u32) (n_dot_l * d * 0xFF), 0xFF);

        // max light
        if (l == 0xFF) {
            break;
        }
    }
    return l;
}

static void calculate_wall_lightmap(
    struct Wall *wall,
    const struct Light *lights,
    usize n_lights) {

    const struct Sector *sector = &state.sectors.arr[wall->sector];

    wall->lightmap.size = (v2i) {
        (int) max(ceilf(wall->len * LIGHTMAP_SCALE), 1),
        (int) max(ceilf(wall->height * LIGHTMAP_SCALE), 1)
    };

    wall->lightmap.data =
        malloc(wall->lightmap.size.x * wall->lightmap.size.y * sizeof(u8));

    for (int j = 0; j < wall->lightmap.size.y; j++) {
        for (int i = 0; i < wall->lightmap.size.x; i++) {
            const f32
                u = clamp(
                    (i + 0.5f) / (wall->len * LIGHTMAP_SCALE), 0.0f, 1.0f),
                v = clamp(
                    (j + 0.5f) / (wall->height * LIGHTMAP_SCALE), 0.0f, 1.0f);

            const v2 p =
                sector_clamp_point(
                    sector,
                    (v2) {
                        wall->a.x + u * wall->d.x,
                        wall->a.y + u * wall->d.y,
                    });

            wall->lightmap.data[j * wall->lightmap.size.x + i] =
                calculate_point_light(
                    lights,
                    n_lights,
                    wall->sector,
                    (v3) {
                        p.x,
                        p.y,
                        state.sectors.arr[wall->sector].floor.z
                            + (v * wall->height)
                    },
                    (v3) {
                        wall->normal.x, wall->normal.y, 0.0f
                    });
        }
    }
}

static void calculate_plane_lightmap(
    const struct Sector *sector,
    struct Plane *plane,
    const struct Light *lights,
    usize n_lights) {

    plane->lightmap.size = (v2i) {
        (int) max(ceilf(sector->span.x * LIGHTMAP_SCALE), 1),
        (int) max(ceilf(sector->span.y * LIGHTMAP_SCALE), 1)
    };

    plane->lightmap.data =
        malloc(plane->lightmap.size.x * plane->lightmap.size.y * sizeof(u8));

    // TODO: a bunch of these points probably aren't in the sector, skip those.
    for (int j = 0; j < plane->lightmap.size.y; j++) {
        for (int i = 0; i < plane->lightmap.size.x; i++) {
            const f32
                u = clamp(
                    (i + 0.5f) / (sector->span.x * LIGHTMAP_SCALE), 0.0f, 1.0f),
                v = clamp(
                    (j + 0.5f) / (sector->span.y * LIGHTMAP_SCALE), 0.0f, 1.0f);

            const v2 p =
                sector_clamp_point(
                    sector,
                    (v2) {
                        sector->min.x + u * sector->span.x,
                        sector->min.y + v * sector->span.y
                    });

            plane->lightmap.data[j * plane->lightmap.size.x + i] =
                calculate_point_light(
                    lights,
                    n_lights,
                    sector->id,
                    (v3) { p.x, p.y, plane->z },
                    plane->type == PLANE_TYPE_FLOOR ?
                        (v3) { 0, 0, 1 }
                        : (v3) { 0, 0, -1 });
        }
    }
}

// lightmap simple sampling at (u, v)
static u8 lightmap_sample(
    u8 *data,
    v2i size,
    f32 u,
    f32 v) {
    u = clamp(u, 0.0f, 0.999f);
    v = clamp(v, 0.0f, 0.999f);
    return data[(((int) (v * size.y)) * size.x) + ((int) (u * size.x))];
}

// lightmap bilinear sampling which uses weighted mean of nearest 4 samples
static u8 lightmap_sample_bilinear(
    u8 *data,
    v2i size,
    f32 u,
    f32 v) {

    const f32
        s_u = 1.0f / size.x,
        s_v = 1.0f / size.y,
        u_r = u - fmodf(u, s_u),
        v_r = v - fmodf(v, s_v),
        inv_d = 1.0f / (s_u * s_v),
        u0 = u_r,
        u1 = u_r + s_u,
        v0 = v_r,
        v1 = v_r + s_v,
        w00 = ((u1 - u) * (v1 - v)) * inv_d,
        w01 = ((u1 - u) * (v - v0)) * inv_d,
        w10 = ((u - u0) * (v1 - v)) * inv_d,
        w11 = ((u - u0) * (v - v0)) * inv_d;

    const u8
        f00 = lightmap_sample(data, size, u0, v0),
        f01 = lightmap_sample(data, size, u0, v1),
        f10 = lightmap_sample(data, size, u1, v0),
        f11 = lightmap_sample(data, size, u1, v1);

    return
        (u8) clamp(
            ((w00 * f00) + (w01 * f01) + (w10 * f10) + (w11 * f11)),
            0x00,
            0xFF);
}

static u8 sample_wall_lightmap(
    const struct Wall *wall,
    f32 u,
    f32 v) {
    return
        lightmap_sample_bilinear(
            wall->lightmap.data, wall->lightmap.size, u, v);
}

static u8 sample_plane_lightmap(
    const struct Sector *sector,
    const struct Plane *plane,
    f32 u,
    f32 v) {
    return
        lightmap_sample_bilinear(
            plane->lightmap.data, plane->lightmap.size, u, v);
}
#endif // ifdef LIGHTMAPS

ALWAYS_INLINE u32 texture_sample_px(
    const struct Texture *tex,
    v2i p) {
    return
        tex->data[
            ((p.y + tex->offset.y) * (tex->pitch / 4))
                + (tex->offset.x + p.x)];
}

ALWAYS_INLINE u32 texture_sample_uv(
    const struct Texture *tex,
    v2 uv) {
    return
        texture_sample_px(
            tex,
            (v2i) {
                clamp(uv.x, 0.0, 1.0 - UV_EPSILON) * tex->size.x,
                clamp(uv.y, 0.0, 1.0 - UV_EPSILON) * tex->size.y
            });
}

static void verline(int x, int y0, int y1, u32 color) {
    for (int y = y0; y <= y1; y++) {
        state.pixels[y * SCREEN_WIDTH + x] = color;
    }
}

static void wall_line(
    const struct Wall *wall,
    int x,
    int y0,
    int y1,
    int y_floor,
    int y_ceil,
    int t_y_floor,
    int t_y_ceil,
    f32 u,
    const struct Texture *texture,
    u8 light) {
    const int tx = (int) (u * texture->size.x);

    const usize pitch_4 = texture->pitch / 4;
    const f32 yd = y_ceil - y_floor;

    for (int y = y0; y <= y1; y++) {
        const f32 v = clamp((y - y_floor) / yd, 0.0f, 1.0f - UV_EPSILON);

        // TODO: TEXTURE_NONE should use pink texture
        const u32 color =
            texture->data[
                ((((int) (v * (texture->size.y - 1))) + texture->offset.y)
                    * pitch_4)
                    + (tx + texture->offset.x)];

        u8 l = light;

#ifdef LIGHTMAPS
        const f32 tv =
            clamp(
                (y - t_y_floor) / (t_y_ceil - t_y_floor),
                0.0,
                1.0 - UV_EPSILON);
        l = clamp(l + (u32) sample_wall_lightmap(wall, u, tv), 0, LIGHT_MAX);
#endif // ifdef LIGHTMAPS

        state.pixels[y * SCREEN_WIDTH + x] = abgr_mul(color, l / (f32) LIGHT_MAX);
    }
}

static void floor_ceil_line(
    const struct Sector *sector,
    const struct Plane *plane,
    int x,
    int y0,
    int y1,
    v2 d_w,
    v2 pw,
    f32 z,
    const struct Texture *tex,
    v2 offset_tex,
    v2 scale_tex,
    u8 light) {
    for (int y = y0; y <= y1; y++) {
        const v2
            pc = wall_floor_ceiling_to_camera_space(d_w, pw, y, z),
            pl = camera_pos_to_world(pc);

        const v2i pt = {
            (int) (max(pl.x, 0) * FLOOR_CEILING_PX_PER_UNIT),
            (int) (max(pl.y, 0) * FLOOR_CEILING_PX_PER_UNIT)
        };

        const u32 color =
            texture_sample_px(
                tex,
                (v2i) {
                    ((int) ((pt.x * scale_tex.x) + offset_tex.x))
                        % tex->size.x,
                    ((int) ((pt.y * scale_tex.y) + offset_tex.y))
                        % tex->size.y
                });


        u8 l = light;

#ifdef LIGHTMAPS
        const u8 lp =
            sample_plane_lightmap(
                sector,
                plane,
                (pl.x - sector->min.x) / sector->span.x,
                (pl.y - sector->min.y) / sector->span.y);
        l = clamp(l + (u32) lp, 0, LIGHT_MAX);
#endif // ifdef LIGHTMAPS

        state.pixels[y * SCREEN_WIDTH + x] =
            abgr_mul(color,  l / (f32) LIGHT_MAX);
    }
}

// DDA line drawing
static void render_line(v2i p0, v2i p1, u32 color) {
    f32
        dx = p1.x - p0.x,
        dy = p1.y - p0.y;

    const f32 step = fabsf(dx) > fabsf(dy) ? fabsf(dx) : fabsf(dy);

    dx /= step;
    dy /= step;

    f32 x = p0.x, y = p0.y;
    int i = 1;
    while (i <= step) {
        const int ix = x, iy = y;
        if (ix >= 0
            && iy >= 0
            && ix < SCREEN_WIDTH
            && iy < SCREEN_HEIGHT) {
            state.pixels[iy * SCREEN_WIDTH + ix] = color;
        }

        x += dx;
        y += dy;
        i++;
    }
}

static void render_sprite(
    int sector,
    v2 pos,
    f32 z,
    v2 scale,
    v2 offset,
    const struct Texture *texture) {
    const f32
        w = texture->size.x * SPRITE_PX_PER_UNIT * scale.x,
        h = texture->size.y * SPRITE_PX_PER_UNIT * scale.y * 2.0f,
        half_w = w / 2.0f;

    const v2 pos_c = world_pos_to_camera(pos);

    if (pos_c.y <= ZNEAR) {
        // sprite is entirely behind camera
        return;
    }

    // arrange points such that camera is always right of p0 -> p1
    const v2
        p0 = { pos_c.x - half_w, pos_c.y },
        p1 = { pos_c.x + half_w, pos_c.y };

    // TODO: can maybe optimize atan2 since "y" is fixed for both points
    // convert to screen angle
    const f32
        ap0 = normalize_angle(atan2(p0.y, p0.x) - PI_2),
        ap1 = normalize_angle(atan2(p1.y, p1.x) - PI_2);

    // check if entirely out of range (both ap0, ap1 are outside of HFOV on the
    // same side)
    if ((ap0 > +(HFOV / 2) && ap1 > (HFOV / 2))
        || (ap0 < -(HFOV / 2) && ap1 < -(HFOV / 2))) {
        return;
    }

    // convert to screen x unbounded
    const int
        tx0 = screen_angle_to_unbounded_x(ap0),
        tx1 = screen_angle_to_unbounded_x(ap1);

    // clamp (will get clamped further by portals)
    int
        x0 = clamp(tx0, 0, SCREEN_WIDTH - 1),
        x1 = clamp(tx1, 0, SCREEN_WIDTH - 1);

    if (x0 == x1) { return; }

    const u8 light = state.sectors.arr[sector].light;

    // TODO: exit early if sprite is entirely behind wall

    // clip y according to portals
    u16 y_lo[SCREEN_WIDTH], y_hi[SCREEN_WIDTH];
    for (int x = x0; x <= x1; x++) {
        y_lo[x] = 0;
        y_hi[x] = SCREEN_HEIGHT - 1;
    }

    // don't line clip if in same sector
    if (sector == state.camera.sector) { goto draw; }

    // create y clipping bounds according to line functions for each column
    enum { QUEUE_MAX = 32 };
    struct { int sectors[SECTOR_MAX]; usize i, n; } queue =
        { { sector }, 0, 1 };

    // compute portal (s)ector(x) min (0) and max (1)
    int sx0 = SCREEN_WIDTH - 1, sx1 = 0;

    while (queue.n != 0) {
        const int id = queue.sectors[queue.i];
        queue.i = (queue.i + 1) % QUEUE_MAX;
        queue.n--;

        if (id == state.camera.sector) {
            continue;
        }

        const struct SectorClip *sc = &state.sectclip[id];

        sx0 = min(sx0, sc->x0);
        sx1 = max(sx1, sc->x1);

        for (usize i = 0; i < sc->n; i++) {
            const struct PortalClip *pc =
                &sc->arr[i];

            // only get clipping lines which affect this x-range
            if (pc->x0 > x1
                || pc->x1 < x0) {
                continue;
            }

            const int
                xx0 = max(x0, pc->x0),
                xx1 = min(x1, pc->x1);

            for (int x = xx0; x <= xx1; x++) {
                const int
                    yf = pc->mf * (x - pc->x0) + pc->yf0,
                    yc = pc->mc * (x - pc->x0) + pc->yc0;

                y_lo[x] = max(yf, y_lo[x]);
                y_hi[x] = min(yc, y_hi[x]);
            }

            queue.sectors[(queue.i + queue.n) % QUEUE_MAX] = pc->from;
            queue.n++;
        }
    }

    x0 = max(x0, sx0);
    x1 = min(x1, sx1);

    if (x0 == x1) { return; }
draw:;
    const f32 sy = ifnan((VFOV * SCREEN_HEIGHT) / pos_c.y, 1e10f);

    const int
        y0 = (SCREEN_HEIGHT / 2) + (int) (((z + 0) - state.camera.z) * sy),
        y1 = (SCREEN_HEIGHT / 2) + (int) (((z + h) - state.camera.z) * sy);

    for (int x = x0; x <= x1; x++) {
        // skip column if behind wall
        if (state.wall_dist[x] < pos_c.y) {
            continue;
        }

        const f32 u = clamp(((x - tx0) / (f32) (tx1 - tx0)), 0.0f, 1.0f);

        const int
            yy0 = clamp(y0, y_lo[x], y_hi[x]),
            yy1 = clamp(y1, y_lo[x], y_hi[x]);

        for (int y = yy0; y < yy1; y++) {
            const f32 v = 1.0f - clamp((y - y0) / (f32) (y1 - y0), 0.0f, 1.0f);

            const u32 color = texture_sample_uv(texture, (v2) { u, v });

            if (((color >> 24) & 0xFF) == 0xFF) {
                state.pixels[y * SCREEN_WIDTH + x] =
                    abgr_mul(color, light / (f32) LIGHT_MAX);
            }
        }
    }
}

static void render() {
    memsetf32(state.wall_dist, 0.0, SCREEN_WIDTH);
    memset16(state.y_hi, SCREEN_HEIGHT - 1, SCREEN_WIDTH);
    memset16(state.y_lo, 0, SCREEN_WIDTH);
    memset(state.sectdraw, 0, sizeof(state.sectdraw));
    memset(state.sectclip, 0, sizeof(state.sectclip));
    memset(state.sectport, 0, sizeof(state.sectport));

    // calculate edges of near/far planes (looking down +Y axis)
    const v2
        zdl = rotate(((v2) { 0.0f, 1.0f }), +(HFOV / 2.0)),
        zdr = rotate(((v2) { 0.0f, 1.0f }), -(HFOV / 2.0)),
        znl = (v2) { zdl.x * ZNEAR, zdl.y * ZNEAR },
        znr = (v2) { zdr.x * ZNEAR, zdr.y * ZNEAR },
        zfl = (v2) { zdl.x * ZFAR, zdl.y * ZFAR },
        zfr = (v2) { zdr.x * ZFAR, zdr.y * ZFAR };

    enum { QUEUE_MAX = 64 };
    struct QueueEntry { int id; struct PortalClip clip; };

    struct { struct QueueEntry arr[QUEUE_MAX]; usize i; usize n; } queue = {
        {{ state.camera.sector, { .x0 = 0, .x1 = SCREEN_WIDTH - 1 }}},
        0,
        1
    };

    while (queue.n != 0) {
        // grab front of queue
        struct QueueEntry entry = queue.arr[queue.i];
        queue.i = (queue.i + 1) % QUEUE_MAX;
        queue.n--;

        // we need to skip sector A (in relation to sector B) if:
        // * we are coming from sector B
        // * there has already been an A to B portal
        if (bitmap_get(state.sectport, (entry.id * SECTOR_MAX) + entry.clip.from)) {
            if (state.sleepy) {
                LOG("skipping portal from %d to %d (at %d %d)", entry.clip.from, entry.id, entry.clip.x0, entry.clip.x1);
            }
            continue;
        }

        // mark this sector as portaled to from the source
        bitmap_set(state.sectport, (entry.clip.from * SECTOR_MAX) + entry.id);

        // mark this sector as "drawn"
        bitmap_set(state.sectdraw, entry.id);

        // add to clipping data
        TYPEOF(&state.sectclip[0]) entry_sectclip = &state.sectclip[entry.id];

        if (entry_sectclip->n == PORTAL_MAX) {
            WARN("out of clip data space for sector %d", entry.id);
        } else {
            entry_sectclip->arr[entry_sectclip->n++] = entry.clip;
        }

        const struct Sector *sector = &state.sectors.arr[entry.id];

        const struct Texture
            *tex_floor = &state.textures[sector->floor.texture],
            *tex_ceil = &state.textures[sector->ceil.texture];

        for (usize i = 0; i < sector->n_walls; i++) {
            const struct Wall *wall =
                &state.walls.arr[sector->first_wall + i];

            const struct Material *material =
                &state.materials.arr[wall->material];

            // translate relative to player and rotate points around player's view
            const v2
                op0 = world_pos_to_camera(AS_V2(wall->a)),
                op1 = world_pos_to_camera(AS_V2(wall->b));

            // wall clipped pos
            v2 cp0 = op0, cp1 = op1;

            // both are negative -> wall is entirely behind player
            if (cp0.y <= 0 && cp1.y <= 0) {
                if(state.sleepy) { LOG("skipping wall %d, entirely behind camera", wall->index); }
                continue;
            }

            // angle-clip against view frustum
            //
            // IMPORTANT
            // screen angles are -HFOV/2 on the RIGHT and +HFOV/2 on the LEFT
            // (the opposite of increasing x-coordinates)
            f32
                ap0 = normalize_angle(atan2(cp0.y, cp0.x) - PI_2),
                ap1 = normalize_angle(atan2(cp1.y, cp1.x) - PI_2);

            // clip against view frustum if both angles are not clearly within
            // HFOV
            if (cp0.y < ZNEAR
                || cp1.y < ZNEAR
                || ap0 > (HFOV / 2.0f)
                || ap1 < -(HFOV / 2.0f)) {
                const v2
                    il = intersect_segs(cp0, cp1, znl, zfl),
                    ir = intersect_segs(cp0, cp1, znr, zfr);

                // recompute angles if points change
                if (!isnan(il.x)) {
                    cp0 = il;
                    ap0 = normalize_angle(atan2(cp0.y, cp0.x) - PI_2);
                }

                if (!isnan(ir.x)) {
                    cp1 = ir;
                    ap1 = normalize_angle(atan2(cp1.y, cp1.x) - PI_2);
                }
            }

            // check if B's angle is greater than A's angle, it ought to be less
            // to be drawn
            if (ap1 > ap0) {
                if(state.sleepy) { LOG("skipping wall %d, angles wrong side", wall->index); }
                continue;
            }

            if ((ap0 > (HFOV / 2.0f) && ap1 > (HFOV / 2.0f))
                || (ap0 < -(HFOV / 2.0f) && ap1 < -(HFOV / 2.0f))) {
                if(state.sleepy) { LOG("skipping wall %d, angles outside of FOV", wall->index); }
                continue;
            }

            // "true" xs before portal clamping
            const int
                tx0 = screen_angle_to_x(ap0),
                tx1 = screen_angle_to_x(ap1);

            // bounds check against portal window
            if (tx0 > entry.clip.x1) { if(state.sleepy) { LOG("skipping wall %d, x0 > x1", wall->index); } continue; }
            if (tx1 < entry.clip.x0) { if(state.sleepy) { LOG("skipping wall %d, x0 > x1", wall->index); } continue; }

            const struct Texture *tex_wall =
                &state.textures[material->texture_id];

            const int
                x0 = clamp(tx0, entry.clip.x0, entry.clip.x1),
                x1 = clamp(tx1, entry.clip.x0, entry.clip.x1);

            if (state.sleepy) {
                LOG("drawing wall %d from (%f %f) %d to %d (%d %d clipped by %d %d)",
                    wall->index,
                    ap0, ap1,
                    x0, x1,
                    tx0, tx1,
                    entry.clip.x0, entry.clip.x1);
            }

            const f32
                z_floor = sector->floor.z,
                z_ceil = sector->ceil.z,
                nz_floor =
                    wall->portal ? state.sectors.arr[wall->portal].floor.z : 0,
                nz_ceil =
                    wall->portal ? state.sectors.arr[wall->portal].ceil.z : 0;

            const f32
                sy0 = ifnan((VFOV * SCREEN_HEIGHT) / cp0.y, 1e10f),
                sy1 = ifnan((VFOV * SCREEN_HEIGHT) / cp1.y, 1e10f);

            const int
                yf0  = (SCREEN_HEIGHT / 2) + (int) (( z_floor - state.camera.z) * sy0),
                yc0  = (SCREEN_HEIGHT / 2) + (int) (( z_ceil  - state.camera.z) * sy0),
                yf1  = (SCREEN_HEIGHT / 2) + (int) (( z_floor - state.camera.z) * sy1),
                yc1  = (SCREEN_HEIGHT / 2) + (int) (( z_ceil  - state.camera.z) * sy1),
                nyf0 = (SCREEN_HEIGHT / 2) + (int) ((nz_floor - state.camera.z) * sy0),
                nyc0 = (SCREEN_HEIGHT / 2) + (int) ((nz_ceil  - state.camera.z) * sy0),
                nyf1 = (SCREEN_HEIGHT / 2) + (int) ((nz_floor - state.camera.z) * sy1),
                nyc1 = (SCREEN_HEIGHT / 2) + (int) ((nz_ceil  - state.camera.z) * sy1),
                txd  = tx1 - tx0,
                yfd  = yf1 - yf0,
                ycd  = yc1 - yc0,
                nyfd = nyf1 - nyf0,
                nycd = nyc1 - nyc0;

            // compute u boundaries according to wall cutoff
            const f32
                u0 = ((cp0.x - op0.x) / (op1.x - op0.x)),
                u1 = ((cp1.x - op0.x) / (op1.x - op0.x));

            const f32
                u0_z0 = u0 / cp0.y,
                u1_z1 = u1 / cp1.y,
                iz0   = 1.0f / cp0.y,
                iz1   = 1.0f / cp1.y;

            // neighbor ceiling/floor bounds for x0/x1
            int yf_x0, yf_x1, yc_x0, yc_x1;

            // TODO: skip this entire loop if nothing will be rendered by height
            for (int x = x0; x <= x1; x++) {
                // calculate progress along x-axis via tx{0,1} so that walls
                // which are partially cut off due to portal edges still have
                // proper heights
                const f32 xp = ifnan((x - tx0) / (f32) txd, 0);

                // perspective correct texture mapping
                // see en.wikipedia.org/wiki/Texture_mapping
                //
                // whereas affine texture mapping from u0..u1 will look like
                // u_a = (1 - a)u_0 + a(u_1) for 0 <= a <= 1
                //
                // perspective-correct texture mapping is of the form
                // u_a =
                //      ((1 - a)(u_0/z_0) + a(u_1/z_1))
                //          / ((1 - a)(1 / z_0) + a(1 / z_1))
                // for the same 0 <= a <= 1
                //
                // where u_0 and u_1 are our calculated-for-affine u0 and u1,
                // a is xp (progress in x direction) and z_0 and z_1 correspond
                // to the y coordinates of the wall points cp0 and cp1
                // respectively
                //
                // we also precalculate values which are consistent from x0..x1
                //
                // "ux" is this wall's u-coordinate as though its entire length
                // were from u in [0..1] (a wall subsection will have a ux
                // subrange of [0..1])
                const f32 ux =
                    clamp(
                        ((((1.0f - xp) * u0_z0) + (xp * u1_z1))
                            / (((1.0f - xp) * iz0) + (xp * iz1))),
                        0.0f,
                        1.0f - UV_EPSILON);

                // get y coordinates for this x
                const int
                    tyf = (int) (xp * yfd) + yf0,
                    tyc = (int) (xp * ycd) + yc0,
                    yf = clamp(tyf, state.y_lo[x], state.y_hi[x]),
                    yc = clamp(tyc, state.y_lo[x], state.y_hi[x]);

                // compute values needed for ceiling/floor texture, which are
                // based off of the fixed world position
                //
                // compute screen angle (theta) and use it to reconstruct:
                // direction to wall (dtw) by inversing the atan2(cp0/1) done
                // above. this will only get us a direction vector, so to find
                // the true point on the wall we need to intersect the line
                // dtw with the original wall in camera space to get us the
                // original wall point(owp)
                const f32 theta = x_to_screen_angle(x);
                const v2
                    d_w = { fastcos(theta + PI_2), fastsin(theta + PI_2) },
                    pw = intersect_lines(((v2) { 0, 0 }), d_w, cp0, cp1);

                // floor
                if (yf > state.y_lo[x]) {
                    floor_ceil_line(
                        sector,
                        &sector->floor,
                        x,
                        state.y_lo[x],
                        yf,
                        d_w,
                        pw,
                        z_floor,
                        tex_floor,
                        sector->floor.offset,
                        sector->floor.scale,
                        sector->light);
                }

                // ceiling
                if (yc < state.y_hi[x]) {
                    floor_ceil_line(
                        sector,
                        &sector->ceil,
                        x,
                        yc,
                        state.y_hi[x],
                        d_w,
                        pw,
                        z_ceil,
                        tex_ceil,
                        sector->ceil.offset,
                        sector->ceil.scale,
                        sector->light);
                }

                if (wall->portal) {
                    const int
                        tnyf = (int) (xp * nyfd) + nyf0,
                        tnyc = (int) (xp * nycd) + nyc0,
                        nyf = clamp(tnyf, state.y_lo[x], state.y_hi[x]),
                        nyc = clamp(tnyc, state.y_lo[x], state.y_hi[x]);

                    // track portal bounds for sprite masking portal data
                    if (x == x0) { yf_x0 = tyf; yc_x0 = tyc; }
                    if (x == x1) { yf_x1 = tyf; yc_x1 = tyc; }

                    wall_line(
                        wall,
                        x,
                        nyc,
                        yc,
                        tyc,
                        tnyc,
                        tyf,
                        tyc,
                        ux,
                        &state.textures[2],
                        sector->light);

                    wall_line(
                        wall,
                        x,
                        yf,
                        nyf,
                        tnyf,
                        tyf,
                        tyf,
                        tyc,
                        ux,
                        &state.textures[3],
                        sector->light);

                    state.y_hi[x] =
                        clamp(
                            min(min(yc, nyc), state.y_hi[x]),
                            0, SCREEN_HEIGHT - 1);

                    state.y_lo[x] =
                        clamp(
                            max(max(yf, nyf), state.y_lo[x]),
                            0, SCREEN_HEIGHT - 1);
                } else {
                    state.wall_dist[x] = pw.y;

                    wall_line(
                        wall,
                        x,
                        yf,
                        yc,
                        tyf,
                        tyc,
                        tyf,
                        tyc,
                        ux,
                        tex_wall,
                        sector->light);
                }

                if (state.sleepy) {
                    present();
                    SDL_Delay(10);
                }
            }

            if (wall->portal) {
                if(state.sleepy) { LOG("queueing sector %d from %d to %d from wall %d (sector %d)", wall->portal, x0, x1, wall->index, sector->id); }
                ASSERT(queue.n != QUEUE_MAX);
                queue.arr[(queue.i + queue.n) % QUEUE_MAX] =
                    (struct QueueEntry) {
                        .id = wall->portal,
                        .clip = {
                            .from = sector->id,
                            .x0 = x0,
                            .x1 = x1,
                            .yf0 = yf_x0,
                            .yc0 = yc_x0,
                            .mf = (yf_x1 - yf_x0) / (f32) (x1 - x0),
                            .mc = (yc_x1 - yc_x0) / (f32) (x1 - x0)
                        }
                    };
                queue.n++;
            }
        }
    }

    // process clipping data for each sector to get entire x0/x1 range
    for (int i = 0; i < SECTOR_MAX; i++) {
        if (!bitmap_get(state.sectdraw, i)) { continue; }

        struct SectorClip *sc = &state.sectclip[i];
        sc->x0 = SCREEN_WIDTH - 1;
        sc->x1 = 0;

        for (usize j = 0; j < sc->n; j++) {
            sc->x0 = min(sc->x0, sc->arr[j].x0);
            sc->x1 = max(sc->x1, sc->arr[j].x1);

            const struct PortalClip *pc = &sc->arr[j];
            if (i != state.camera.sector) {
                render_line(
                    (v2i) { pc->x0, pc->yf0 },
                    (v2i) { pc->x1, pc->mf * (pc->x1 - pc->x0) + pc->yf0 },
                    0xFFFF00FF);

                render_line(
                    (v2i) { pc->x0, pc->yc0 },
                    (v2i) { pc->x1, pc->mc * (pc->x1 - pc->x0) + pc->yc0 },
                    0xFF00FFFF);

                render_line((v2i) { pc->x0, 0 }, (v2i) { pc->x0, SCREEN_HEIGHT - 1 }, 0xFF0000FF);
                render_line((v2i) { pc->x1, 0 }, (v2i) { pc->x1, SCREEN_HEIGHT - 1 }, 0xFFFF0000);
            }
        }
    }

    render_sprite(
        2,
        (v2) { 7, 6 },
        1.0,
        (v2) { 1, 1 },
        (v2) { 1, 1 },
        &state.textures[TEXTURE_TEST4]);

    state.sleepy = false;
}

static void present() {
    void *px;
    int pitch;
    SDL_LockTexture(state.texture, NULL, &px, &pitch);
    {
        for (usize y = 0; y < SCREEN_HEIGHT; y++) {
            memcpy(
                &((u8*) px)[y * pitch],
                &state.pixels[y * SCREEN_WIDTH],
                SCREEN_WIDTH * 4);
        }
    }
    SDL_UnlockTexture(state.texture);

    SDL_SetRenderTarget(state.renderer, NULL);
    SDL_SetRenderDrawColor(state.renderer, 0, 0, 0, 0xFF);
    SDL_SetRenderDrawBlendMode(state.renderer, SDL_BLENDMODE_NONE);

    SDL_RenderClear(state.renderer);
    SDL_RenderCopyEx(
        state.renderer,
        state.texture,
        NULL,
        NULL,
        0.0,
        NULL,
        SDL_FLIP_VERTICAL);

    SDL_SetTextureBlendMode(state.debug, SDL_BLENDMODE_BLEND);
    SDL_RenderCopy(state.renderer, state.debug, NULL, &((SDL_Rect) { 0, 0, 512, 512 }));
    SDL_RenderPresent(state.renderer);
}

int main(int argc, char *argv[]) {
    make_trigtabs();

    ASSERT(
        !SDL_Init(SDL_INIT_VIDEO),
        "SDL failed to initialize: %s",
        SDL_GetError());

    ASSERT(
        IMG_Init(IMG_INIT_PNG) & IMG_INIT_PNG,
        "failed to load SDL_image: %s",
        IMG_GetError());

    state.window =
        SDL_CreateWindow(
            "raycast",
            SDL_WINDOWPOS_CENTERED_DISPLAY(1),
            SDL_WINDOWPOS_CENTERED_DISPLAY(1),
            1280,
            720,
            0);

    ASSERT(state.window, "failed to create SDL window: %s\n", SDL_GetError());

    state.renderer =
        SDL_CreateRenderer(
            state.window,
            -1,
            SDL_RENDERER_ACCELERATED
            | SDL_RENDERER_PRESENTVSYNC);

    state.texture =
        SDL_CreateTexture(
            state.renderer,
            SDL_PIXELFORMAT_ABGR8888,
            SDL_TEXTUREACCESS_STREAMING,
            SCREEN_WIDTH,
            SCREEN_HEIGHT);
    state.debug =
        SDL_CreateTexture(
            state.renderer,
            SDL_PIXELFORMAT_ABGR8888,
            SDL_TEXTUREACCESS_TARGET,
            128,
            128);

    state.pixels = malloc(SCREEN_WIDTH * SCREEN_HEIGHT * 4);

    SDL_Surface
        *texture_surface_loaded = IMG_Load("res/test.png"),
        *texture_surface =
            SDL_ConvertSurfaceFormat(
                texture_surface_loaded,
                SDL_PIXELFORMAT_ABGR8888,
                0);
    state.texture_pixels = texture_surface->pixels;

    state.textures[TEXTURE_TEST0] = (struct Texture) {
        .id = TEXTURE_TEST0,
        .data = state.texture_pixels,
        .offset = { 0, 0 },
        .size = { 16, 16 },
        .pitch = texture_surface->pitch
    };

    state.textures[TEXTURE_TEST1] = (struct Texture) {
        .id = TEXTURE_TEST1,
        .data = state.texture_pixels,
        .offset = { 16, 0 },
        .size = { 16, 16 },
        .pitch = texture_surface->pitch
    };


    state.textures[TEXTURE_TEST2] = (struct Texture) {
        .id = TEXTURE_TEST2,
        .data = state.texture_pixels,
        .offset = { 32, 0 },
        .size = { 16, 16 },
        .pitch = texture_surface->pitch
    };


    state.textures[TEXTURE_TEST3] = (struct Texture) {
        .id = TEXTURE_TEST3,
        .data = state.texture_pixels,
        .offset = { 48, 0 },
        .size = { 16, 16 },
        .pitch = texture_surface->pitch
    };

    state.textures[TEXTURE_TEST4] = (struct Texture) {
        .id = TEXTURE_TEST4,
        .data = state.texture_pixels,
        .offset = { 64, 0 },
        .size = { 16, 16 },
        .pitch = texture_surface->pitch
    };

    state.camera.pos = (v2) { 3.04, 3 };
    state.camera.angle = 0.0;
    state.camera.sector = 1;

    int ret = 0;
    ASSERT(
        !(ret = load_sectors("res/level.txt")),
        "error while loading sectors: %d",
        ret);
    LOG(
        "loaded %" PRIusize " sectors with %" PRIusize " walls",
        state.sectors.n,
        state.walls.n);

    state.materials.arr[state.materials.n++] = (struct Material) {
        .id = state.materials.n - 1,
        .texture_id = TEXTURE_TEST0,
        .scale = { 1, 1 }
    };

    // calculate lighting data
#ifdef LIGHTMAPS
    for (usize i = SECTOR_FIRST; i < state.sectors.n; i++) {
        struct Sector *sector = &state.sectors.arr[i];
        calculate_plane_lightmap(sector, &sector->floor, g_lights, ARRLEN(g_lights));
        calculate_plane_lightmap(sector, &sector->ceil, g_lights, ARRLEN(g_lights));

        for (usize j = 0; j < sector->n_walls; j++) {
            calculate_wall_lightmap(
                &state.walls.arr[sector->first_wall + j], g_lights, ARRLEN(g_lights));
        }
    }
#endif

    state.camera.z = EYE_Z;

    while (!state.quit) {
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) {
            switch (ev.type) {
                case SDL_QUIT:
                    state.quit = true;
                    break;
                default:
                    break;
            }
        }

        if (state.quit) {
            break;
        }

        const f32 rot_speed = 3.0f * 0.032f, move_speed = 3.0f * 0.032f;

        const u8 *keystate = SDL_GetKeyboardState(NULL);

        if (keystate[SDLK_RIGHT & 0xFFFF]) {
            state.camera.angle -= rot_speed;
        }

        if (keystate[SDLK_LEFT & 0xFFFF]) {
            state.camera.angle += rot_speed;
        }

        state.camera.anglecos = fastcos(state.camera.angle);
        state.camera.anglesin = fastsin(state.camera.angle);

        if (keystate[SDLK_UP & 0xFFFF]) {
            state.camera.pos = (v2) {
                state.camera.pos.x + (move_speed * state.camera.anglecos),
                state.camera.pos.y + (move_speed * state.camera.anglesin),
            };
        }

        if (keystate[SDLK_DOWN & 0xFFFF]) {
            state.camera.pos = (v2) {
                state.camera.pos.x - (move_speed * state.camera.anglecos),
                state.camera.pos.y - (move_speed * state.camera.anglesin),
            };
        }

        if (keystate[SDLK_F1 & 0xFFFF]) {
            state.sleepy = true;
        }

        // update player sector
        {
            // BFS neighbors in a circular queue, player is likely to be in one
            // of the neighboring sectors
            enum { QUEUE_MAX = 64 };
            int
                queue[QUEUE_MAX] = { state.camera.sector },
                i = 0,
                n = 1,
                found = SECTOR_NONE;

            while (n != 0) {
                // get front of queue and advance to next
                const int id = queue[i];
                i = (i + 1) % (QUEUE_MAX);
                n--;

                const struct Sector *sector = &state.sectors.arr[id];

                if (point_in_sector(sector, state.camera.pos)) {
                    found = id;
                    break;
                }

                // check neighbors
                for (usize j = 0; j < sector->n_walls; j++) {
                    const struct Wall *wall =
                        &state.walls.arr[sector->first_wall + j];

                    if (wall->portal) {
                        if (n == QUEUE_MAX) {
                            WARN("out of queue space!");
                            goto done;
                        }

                        queue[(i + n) % QUEUE_MAX] = wall->portal;
                        n++;
                    }
                }
            }


done:
            if (!found) {
                WARN("player is not in a sector!");
                state.camera.sector = 1;
            } else {
                state.camera.sector = found;
            }
        }

        state.camera.z =
            state.sectors.arr[state.camera.sector].floor.z + EYE_Z;

        memset(state.pixels, 0, SCREEN_WIDTH * SCREEN_HEIGHT * 4);
        render();
        if (!state.sleepy) { present(); }
    }

    SDL_FreeSurface(texture_surface_loaded);
    SDL_FreeSurface(texture_surface);
    SDL_DestroyTexture(state.debug);
    SDL_DestroyTexture(state.texture);
    SDL_DestroyRenderer(state.renderer);
    SDL_DestroyWindow(state.window);
    return 0;
}
