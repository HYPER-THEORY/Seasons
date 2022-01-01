#ifndef WINDOW_HPP
#define WINDOW_HPP

#include <cmath>
#include <string>
#include <vector>
#include <SDL2/SDL.h>

class window {
public:
	static bool isopen;
	static uint32_t delay;
	static std::vector<uint32_t> canvas;
	static bool keydown[128];
	static bool keypressed[128];
	static bool keyreleased[128];
	static bool mouselock;
	static int mousex;
	static int mousey;
	
	static void init(std::string t, int w, int h, int x = SDL_WINDOWPOS_CENTERED, int y = SDL_WINDOWPOS_CENTERED) {
		width = w;
		height = h;
		SDL_Init(SDL_INIT_VIDEO);
		sdlwindow = SDL_CreateWindow(t.c_str(), x, y, w, h, SDL_WINDOW_SHOWN);
		surface = SDL_GetWindowSurface(sdlwindow);
		canvas = std::vector<uint32_t>(w * h);
		std::fill(keydown, keydown + 128, false);
		std::fill(keypressed, keypressed + 128, false);
		std::fill(keyreleased, keyreleased + 128, false);
	}
	
	static void close() {
		isopen = false;
		SDL_DestroyWindow(sdlwindow);
		SDL_Quit();
	}
	
	static uint32_t update() {
		std::fill(keypressed, keypressed + 128, false);
		std::fill(keyreleased, keyreleased + 128, false);
		SDL_Event event;
		while (SDL_PollEvent(&event)) {
			int32_t keycode = event.key.keysym.sym;
			switch (event.type) {
				case SDL_QUIT:
					close();
					return 0;
				case SDL_KEYDOWN:
					if (keycode < 128) {
						keypressed[keycode] = !keydown[keycode];
						keydown[keycode] = true;
					}
					break;
				case SDL_KEYUP:
					if (keycode < 128) {
						keyreleased[keycode] = true;
						keydown[keycode] = false;
					}
					break;
				case SDL_MOUSEBUTTONDOWN:
					if (event.button.button == SDL_BUTTON_LEFT) {
						keypressed[1] = true;
						keydown[1] = true;
					} else if (event.button.button == SDL_BUTTON_RIGHT) {
						keypressed[2] = true;
						keydown[2] = true;
					}
					break;
				case SDL_MOUSEBUTTONUP:
					if (event.button.button == SDL_BUTTON_LEFT) {
						keypressed[1] = true;
						keydown[1] = false;
					} else if (event.button.button == SDL_BUTTON_RIGHT) {
						keypressed[2] = true;
						keydown[2] = false;
					}
					break;
				case SDL_MOUSEMOTION:
					mousex = event.motion.x;
					mousey = event.motion.y;
					break;
			}
		}
		if (mouselock && SDL_GetKeyboardFocus()) SDL_WarpMouseInWindow(sdlwindow, width / 2, height / 2);
		copy(canvas.begin(), canvas.end(), static_cast<uint32_t*>(surface->pixels));
		SDL_UpdateWindowSurface(sdlwindow);
		uint32_t deltatime = SDL_GetTicks() - time;
		if (deltatime <= delay) SDL_Delay(delay - deltatime);
		time = SDL_GetTicks();
		return std::max(delay, deltatime);
	}
	
	static void mousemove(int x, int y) {
		SDL_WarpMouseInWindow(sdlwindow, x, y);
		mousex = x;
		mousey = y;
	}
	
	static void showcursor() {
		SDL_ShowCursor(SDL_ENABLE);
	}
	
	static void hidecursor() {
		SDL_ShowCursor(SDL_DISABLE);
	}
	
private:
	static int width;
	static int height;
	static uint32_t time;
	static SDL_Window* sdlwindow;
	static SDL_Surface* surface;
};

inline int window::width = 0;
inline int window::height = 0;
inline uint32_t window::time = 0;
inline SDL_Window* window::sdlwindow = nullptr;
inline SDL_Surface* window::surface = nullptr;
inline bool window::isopen = true;
inline std::uint32_t window::delay = 0;
inline std::vector<uint32_t> window::canvas;
inline bool window::keydown[128];
inline bool window::keypressed[128];
inline bool window::keyreleased[128];
inline bool window::mouselock = false;
inline int window::mousex = 0;
inline int window::mousey = 0;

#endif
