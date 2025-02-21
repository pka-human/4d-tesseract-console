#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#include <conio.h>
#else
#include <sys/ioctl.h>
#include <termios.h>
#include <fcntl.h>
#endif

typedef struct {
    uint8_t x, y;
} Vector2;

typedef struct {
    int8_t x, y, z;
} Vector3;

typedef struct {
    int8_t x, y, z, w;
} Vector4;

typedef struct {
    Vector4 a, b;
} drawing;

drawing *drawings = NULL;
drawing *drawings_buffer = NULL;
size_t drawings_size = 0;

uint8_t* screen = NULL;
uint8_t screen_x;
uint8_t screen_y;

float PIXEL_ASPECT = 1.0f;

unsigned previous_rows;
unsigned previous_cols;

long long get_microseconds() {
#ifdef _WIN32
    LARGE_INTEGER freq, time;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&time);
    return (time.QuadPart * 1000000LL) / freq.QuadPart;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000000LL + ts.tv_nsec / 1000;
#endif
}

void get_terminal_size(unsigned *rows, unsigned *cols) {
#ifdef _WIN32
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
    *rows = csbi.srWindow.Bottom - csbi.srWindow.Top + 1;
    *cols = csbi.srWindow.Right - csbi.srWindow.Left + 1;
#else
    struct winsize w;
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) == -1) {
        *rows = 0;
        *cols = 0;
    } else {
        *rows = w.ws_row;
        *cols = w.ws_col;
    }
#endif
}

void clear_terminal() {
#ifdef _WIN32
    HANDLE hOut;
    COORD Position;
    hOut = GetStdHandle(STD_OUTPUT_HANDLE);
    Position.X = 0;
    Position.Y = 0;
    SetConsoleCursorPosition(hOut, Position);
#else
    printf("\033[2J");
    printf("\033[H");
#endif
}

float get_char_aspect_ratio() {
    unsigned height = 5;
    unsigned width = height;
    float aspect_ratio = 1.0;
    int ch;
    bool changed = true;

    #ifdef _WIN32
        HANDLE hStdin = GetStdHandle(STD_INPUT_HANDLE);
        DWORD mode;
        GetConsoleMode(hStdin, &mode);
        SetConsoleMode(hStdin, mode & ~(ENABLE_LINE_INPUT | ENABLE_ECHO_INPUT));
    #else
        struct termios oldt, newt;
        tcgetattr(STDIN_FILENO, &oldt);
        newt = oldt;
        newt.c_lflag &= ~(ICANON | ECHO);
        tcsetattr(STDIN_FILENO, TCSANOW, &newt);
        fcntl(STDIN_FILENO, F_SETFL, fcntl(STDIN_FILENO, F_GETFL) | O_NONBLOCK);
    #endif

    while (1) {
        if (changed) {
            #ifdef _WIN32
                system("cls");
            #else
                clear_terminal();
            #endif
            printf("Use left/right arrow keys to adjust\nthe width until it looks like a square.\nPress Enter when it's right.\n\n");

            for (unsigned i = 0; i < height; ++i) {
                for (unsigned j = 0; j < width; ++j) {
                    putchar((i == 0 || i == height-1 || j == 0 || j == width-1) ? '@' : '.');
                }
                putchar('\n');
            }
            changed = false;
        }

        #ifdef _WIN32
            if (_kbhit()) {
                ch = _getch();
                if (ch == 0xE0 || ch == 0x00) {
                    ch = _getch();
                    switch(ch) {
                        case 'K': if(width > height) { --width; changed = true; } break;
                        case 'M': if(width < 15) { ++width; changed = true; } break;
                    }
                }
                else if (ch == '\r') break;
            }
        #else
            if ((ch = getchar()) != EOF) {
                if (ch == '\n' || ch == '\r') break;
                if (ch == '\033') {
                    getchar();
                    switch(getchar()) {
                        case 'D': if(width > height) { --width; changed = true; } break;
                        case 'C': if(width < 15) { ++width; changed = true; } break;
                    }
                }
            }
        #endif
        usleep(10000);
    }

    #ifdef _WIN32
        SetConsoleMode(hStdin, mode);
        system("cls");
    #else
        tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
    #endif
    return (float)height / width;
}

void set_bit(uint8_t x, uint8_t y, bool value) {
    uint16_t pos = y * screen_x + x;
    uint16_t byte_idx = pos / 8;
    uint8_t bit_offset = pos % 8;

    if (value) {
        screen[byte_idx] |= (1 << (7 - bit_offset));
    } else {
        screen[byte_idx] &= ~(1 << (7 - bit_offset));
    }
}

bool get_bit(uint8_t x, uint8_t y) {
    uint16_t pos = y * screen_x + x;
    uint16_t byte_idx = pos / 8;
    uint8_t bit_offset = pos % 8;
    return (screen[byte_idx] >> (7 - bit_offset)) & 1;
}

void init_screen() {
    size_t bytes_needed = (screen_x * screen_y + 7) / 8;
    screen = (uint8_t*)malloc(bytes_needed);
    if (screen == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for screen!\n");
        exit(1);
    }
}

void reinit_screen() {
    free(screen);
    screen = NULL;
    init_screen();
}

void calculate_screen_size(unsigned terminal_rows, unsigned terminal_cols) {
  if (terminal_rows == 0 || terminal_cols == 0) {
    screen_x = 0;
    screen_y = 0;
    return;
  }

  unsigned max_y = terminal_rows < 4 ? 0 : fmin(terminal_rows - 3, 255.0);
  unsigned max_x = fmin(terminal_cols - 1, 255.0);

  unsigned potential_x = round(max_y / PIXEL_ASPECT);
  unsigned potential_y = round(max_x * PIXEL_ASPECT);

  if (potential_x <= max_x && potential_y <= max_y)
   {
    screen_x = (uint8_t) potential_x;
    screen_y = (uint8_t) potential_y;
   } 
  else if (potential_x <= max_x)
   {
      screen_x = (uint8_t) potential_x;
      screen_y = (uint8_t) fmin(round(screen_x * PIXEL_ASPECT), max_y);
   }
  else if(potential_y <= max_y)
  {
    screen_y = (uint8_t) potential_y;
    screen_x = (uint8_t) fmin(round(screen_y / PIXEL_ASPECT), max_x);
  }
  else 
   {
        if (potential_x > max_x) {
            screen_x = (uint8_t)max_x;
            screen_y = (uint8_t)fmin(round(screen_x * PIXEL_ASPECT), max_y);
        } else {
            screen_y = (uint8_t)max_y;
            screen_x = (uint8_t)fmin(round(screen_y / PIXEL_ASPECT), max_x);
        }
    }
}

bool update_screen_size() {
    unsigned terminal_rows, terminal_cols;

    get_terminal_size(&terminal_rows, &terminal_cols);

    if (terminal_rows == previous_rows && terminal_cols == previous_cols) {
        return false;
    }

    previous_rows = terminal_rows;
    previous_cols = terminal_cols;

    calculate_screen_size(terminal_rows, terminal_cols);

    return true;
}

void free_all() {
    free(drawings);
    free(drawings_buffer);
    drawings = NULL;
    drawings_buffer = NULL;
    drawings_size = 0;
}

void clear_screen() {
    size_t bytes = (screen_x * screen_y + 7) / 8;
    memset(screen, 0, bytes);
}

void reallocate_drawings_buffer() {
    drawings_buffer = (drawing*) realloc(drawings_buffer, (++drawings_size) * sizeof(drawing));
    if (drawings_buffer == NULL) {
        fprintf(stderr, "Error: Failed to reallocate memory for drawings_buffer!\n");
        exit(1);
    }
}

void allocate_drawings() {
    if (drawings != NULL) {
        return;
    }
    drawings = (drawing*) malloc(drawings_size * sizeof(drawing));
    if (drawings == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for drawings!\n");
        exit(1);
    }
}

Vector2 project3d2d(bool is_perspective, Vector3 point, float fov_degrees, float zoom) {
    Vector2 result;
    float normalized_x, normalized_y;
    float scale;

    normalized_x = ((float)point.x / 100.0f);
    normalized_y = ((float)point.y / 100.0f);

    if (!is_perspective) {
        scale = zoom;
    } else {
        float fov_radians = fov_degrees * M_PI / 180.0f;
        float z_f = (float)point.z;
        float dist = 1.0f/tanf(fov_radians / 2.0f) * 100.0f;

        if (z_f + dist <= 0.0f) {
            scale = dist / (dist - z_f) * zoom;
        } else {
            scale = dist / (z_f + dist) * zoom;
        }
    }

    int result_x = (int)round((normalized_x * scale + 1.0f) / 2.0f * (float)(screen_x - 1));
    if (result_x >= screen_x) {
        result_x = screen_x - 1;
    } else if (result_x < 0) {
        result_x = 0;
    }

    int result_y = (int)round((normalized_y * scale + 1.0f) / 2.0f * (float)(screen_y - 1));
    if (result_y >= screen_y) {
        result_y = screen_y - 1;
    } else if (result_y < 0) {
        result_y = 0;
    }

    result.x = result_x;
    result.y = result_y;

    return result;
}

Vector3 project4d3d(bool is_perspective, Vector4 point, float fov_degrees, float zoom) {
    Vector3 result;
    float scale;

    if (!is_perspective) {
        scale = zoom;
    } else {
        float fov_radians = fov_degrees * M_PI / 180.0f;
        float w_f = (float)point.w;
        float dist = 1.0f / tanf(fov_radians / 2.0f) * 100.0f;

        if (w_f + dist <= 0.0f) {
            scale = dist / (dist - w_f) * zoom;
        } else {
            scale = dist / (w_f + dist) * zoom;
        }
    }

    result.x = (float)point.x * scale;
    result.y = (float)point.y * scale;
    result.z = (float)point.z * scale;

    return result;
}

void line(Vector4 point_a, Vector4 point_b) {
    if (point_a.x >= -100 && point_a.x <= 100 &&
        point_a.y >= -100 && point_a.y <= 100 &&
        point_a.z >= -100 && point_a.z <= 100 &&
        point_a.w >= -100 && point_a.w <= 100 &&
        point_b.x >= -100 && point_b.x <= 100 &&
        point_b.y >= -100 && point_b.y <= 100 &&
        point_b.z >= -100 && point_b.z <= 100 &&
        point_b.w >= -100 && point_b.w <= 100) {

        reallocate_drawings_buffer();

        drawings_buffer[drawings_size - 1] = (drawing) {
            .a = point_a,
            .b = point_b
        };
    } else {
        fprintf(stderr, "Error: Trying to draw a line or part of a line out of range!\n");
        exit(1);
    }
}

Vector4 rotateXY(Vector4 v, float angle) {
    float cosA = cosf(angle);
    float sinA = sinf(angle);
    return (Vector4){ 
        v.x * cosA - v.y * sinA, 
        v.x * sinA + v.y * cosA, 
        v.z, 
        v.w 
    };
}

Vector4 rotateXZ(Vector4 v, float angle) {
    float cosA = cosf(angle);
    float sinA = sinf(angle);
    return (Vector4){ 
        v.x * cosA - v.z * sinA, 
        v.y, 
        v.x * sinA + v.z * cosA, 
        v.w 
    };
}

Vector4 rotateXW(Vector4 v, float angle) {
    float cosA = cosf(angle);
    float sinA = sinf(angle);
    return (Vector4){ 
        v.x * cosA - v.w * sinA, 
        v.y, 
        v.z, 
        v.x * sinA + v.w * cosA 
    };
}

Vector4 rotateYZ(Vector4 v, float angle) {
    float cosA = cosf(angle);
    float sinA = sinf(angle);
    return (Vector4){ 
        v.x, 
        v.y * cosA - v.z * sinA, 
        v.y * sinA + v.z * cosA, 
        v.w 
    };
}

Vector4 rotateYW(Vector4 v, float angle) {
    float cosA = cosf(angle);
    float sinA = sinf(angle);
    return (Vector4){ 
        v.x, 
        v.y * cosA - v.w * sinA, 
        v.z, 
        v.y * sinA + v.w * cosA 
    };
}

Vector4 rotateZW(Vector4 v, float angle) {
    float cosA = cosf(angle);
    float sinA = sinf(angle);
    return (Vector4){ 
        v.x, 
        v.y, 
        v.z * cosA - v.w * sinA, 
        v.z * sinA + v.w * cosA 
    };
}

void rotate_world_XY(float theta) {
    for (size_t i = 0; i < drawings_size; ++i) {
        drawings[i].a = rotateXY(drawings[i].a, theta);
        drawings[i].b = rotateXY(drawings[i].b, theta);
    }
}

void rotate_world_XZ(float theta) {
    for (size_t i = 0; i < drawings_size; ++i) {
        drawings[i].a = rotateXZ(drawings[i].a, theta);
        drawings[i].b = rotateXZ(drawings[i].b, theta);
    }
}

void rotate_world_XW(float theta) {
    for (size_t i = 0; i < drawings_size; ++i) {
        drawings[i].a = rotateXW(drawings[i].a, theta);
        drawings[i].b = rotateXW(drawings[i].b, theta);
    }
}

void rotate_world_YZ(float theta) {
    for (size_t i = 0; i < drawings_size; ++i) {
        drawings[i].a = rotateYZ(drawings[i].a, theta);
        drawings[i].b = rotateYZ(drawings[i].b, theta);
    }
}

void rotate_world_YW(float theta) {
    for (size_t i = 0; i < drawings_size; ++i) {
        drawings[i].a = rotateYW(drawings[i].a, theta);
        drawings[i].b = rotateYW(drawings[i].b, theta);
    }
}

void rotate_world_ZW(float theta) {
    for (size_t i = 0; i < drawings_size; ++i) {
        drawings[i].a = rotateZW(drawings[i].a, theta);
        drawings[i].b = rotateZW(drawings[i].b, theta);
    }
}

void draw_line2d(Vector2 point_a, Vector2 point_b) {
    int dx = abs(point_b.x - point_a.x);
    int dy = abs(point_b.y - point_a.y);
    int sx = (point_a.x < point_b.x) ? 1 : -1;
    int sy = (point_a.y < point_b.y) ? 1 : -1;
    int err = dx - dy;
    int e2;

    while (true) {
        if (point_a.x >= 0 && point_a.x < screen_x && point_a.y >= 0 && point_a.y < screen_y) {
            set_bit(point_a.x, point_a.y, true);
        }

        if (point_a.x == point_b.x && point_a.y == point_b.y)
            break;

        e2 = 2 * err;
        if (e2 > -dy) {
            err -= dy;
            point_a.x += sx;
        }
        if (e2 < dx) {
            err += dx;
            point_a.y += sy;
        }

        if (e2 <= -dy && e2 >= dx) {
            if (dx > dy) {
                point_a.x += sx;
            } else {
                point_a.y += sy;
            }
        }
    }
}

void draw(bool perspective, float fov_degrees, float zoom) {
    clear_screen();

    if (screen_x == 0 || screen_y == 0) {
        return;
    }

    if (previous_cols <= 38 || previous_rows <= 15) {
        printf("Too small :(\n");
        return;
    }

    printf("4D Tesseract in console. (Ctrl + C to quit)\nWritten in C by pka_human, 2025.\n");

    for (size_t i = 0; i < drawings_size; ++i) {
        drawing d = drawings[i];
        Vector3 projected_a = project4d3d(perspective, d.a, fov_degrees, zoom);
        Vector3 projected_b = project4d3d(perspective, d.b, fov_degrees, zoom);
        Vector2 draw_point_a = project3d2d(perspective, projected_a, fov_degrees, zoom);
        Vector2 draw_point_b = project3d2d(perspective, projected_b, fov_degrees, zoom);
        draw_line2d(draw_point_a, draw_point_b);
    }

    char line_buf[screen_x + 3]; 

    for (uint8_t yp = 0; yp < screen_y; ++yp) {
        size_t pos = 0;
        for (uint8_t xp = 0; xp < screen_x; ++xp) {
            line_buf[pos++] = get_bit(xp, yp) ? '@' : ' ';
        }
        line_buf[pos++] = '|';
        line_buf[pos++] = '\n';
        line_buf[pos] = '\0';
        printf("%s", line_buf);
    }
}


void tesseract(const int8_t s) {
    line((Vector4){-s, -s, -s, -s}, (Vector4){ s, -s, -s, -s});
    line((Vector4){ s, -s, -s, -s}, (Vector4){ s,  s, -s, -s});
    line((Vector4){ s,  s, -s, -s}, (Vector4){-s,  s, -s, -s});
    line((Vector4){-s,  s, -s, -s}, (Vector4){-s, -s, -s, -s});

    line((Vector4){-s, -s,  s, -s}, (Vector4){ s, -s,  s, -s});
    line((Vector4){ s, -s,  s, -s}, (Vector4){ s,  s,  s, -s});
    line((Vector4){ s,  s,  s, -s}, (Vector4){-s,  s,  s, -s});
    line((Vector4){-s,  s,  s, -s}, (Vector4){-s, -s,  s, -s});

    line((Vector4){-s, -s, -s, -s}, (Vector4){-s, -s,  s, -s});
    line((Vector4){ s, -s, -s, -s}, (Vector4){ s, -s,  s, -s});
    line((Vector4){ s,  s, -s, -s}, (Vector4){ s,  s,  s, -s});
    line((Vector4){-s,  s, -s, -s}, (Vector4){-s,  s,  s, -s});

    line((Vector4){-s, -s, -s,  s}, (Vector4){ s, -s, -s,  s});
    line((Vector4){ s, -s, -s,  s}, (Vector4){ s,  s, -s,  s});
    line((Vector4){ s,  s, -s,  s}, (Vector4){-s,  s, -s,  s});
    line((Vector4){-s,  s, -s,  s}, (Vector4){-s, -s, -s,  s});

    line((Vector4){-s, -s,  s,  s}, (Vector4){ s, -s,  s,  s});
    line((Vector4){ s, -s,  s,  s}, (Vector4){ s,  s,  s,  s});
    line((Vector4){ s,  s,  s,  s}, (Vector4){-s,  s,  s,  s});
    line((Vector4){-s,  s,  s,  s}, (Vector4){-s, -s,  s,  s});

    line((Vector4){-s, -s, -s,  s}, (Vector4){-s, -s,  s,  s});
    line((Vector4){ s, -s, -s,  s}, (Vector4){ s, -s,  s,  s});
    line((Vector4){ s,  s, -s,  s}, (Vector4){ s,  s,  s,  s});
    line((Vector4){-s,  s, -s,  s}, (Vector4){-s,  s,  s,  s});

    line((Vector4){-s, -s, -s, -s}, (Vector4){-s, -s, -s,  s});
    line((Vector4){ s, -s, -s, -s}, (Vector4){ s, -s, -s,  s});
    line((Vector4){ s,  s, -s, -s}, (Vector4){ s,  s, -s,  s});
    line((Vector4){-s,  s, -s, -s}, (Vector4){-s,  s, -s,  s});

    line((Vector4){-s, -s,  s, -s}, (Vector4){-s, -s,  s,  s});
    line((Vector4){ s, -s,  s, -s}, (Vector4){ s, -s,  s,  s});
    line((Vector4){ s,  s,  s, -s}, (Vector4){ s,  s,  s,  s});
    line((Vector4){-s,  s,  s, -s}, (Vector4){-s,  s,  s,  s});
}

int main() {
    PIXEL_ASPECT = get_char_aspect_ratio();

    update_screen_size();
    init_screen();

    long long previous_time;

    float rotation = 0;

    tesseract(50); // tesseract, size 50.

    allocate_drawings();

    bool perspective = true;
    float fov_degrees = 60.0f;
    float zoom = 0.8f;

    while (1) {
        previous_time = get_microseconds();

        memcpy(drawings, drawings_buffer, drawings_size * sizeof(drawing));

        rotate_world_XZ(rotation);
        rotate_world_YW(rotation);
        
        if (update_screen_size()) {
            reinit_screen();
            #ifdef _WIN32
                system("cls");
            #endif
        }

        clear_terminal();
        draw(perspective, fov_degrees, zoom);

        rotation += 0.008f;

        long long sleep = 10000 - (get_microseconds() - previous_time);
        if (sleep > 0) {
            usleep(sleep);
        }
    }

    free_all();
    return 0;
}

