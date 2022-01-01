#ifndef GAME_HPP
#define GAME_HPP

#include <string>
#include <vector>
#include "window.hpp"

void load();
void update(double dt);

inline int width = 960;
inline int height = 540;
inline int scale = 1;
inline std::string title;
inline int fps = 30;
inline bool hidecursor = false;
inline bool firstperson = false;
inline std::vector<uint32_t> framebuffer;

inline void run() {
	load();
	// initialize window
	window::init(title, width * scale, height * scale);
	window::delay = 1000 / fps;
	if (hidecursor) window::hidecursor();
	if (firstperson) {
		window::mouselock = true;
		window::mousemove(width * scale / 2, height * scale / 2);
	}
	// frame buffer
	framebuffer = std::vector<uint32_t>(width * height);
	// delta time
	uint32_t deltatime = 0;
	// main loop time
	while (true) {
		update(deltatime * .001);
		// resize display size
		for (int x = width * scale; x --> 0;)
			for (int y = height * scale; y --> 0;)
				window::canvas[x + y * width * scale] = framebuffer[x / scale + y / scale * width];
		// update window
		deltatime = window::update();
		if (window::keydown[SDLK_ESCAPE]) window::close();
		if (!window::isopen) break;
		printf("FPS: %d\n", 1000 / deltatime);
	}
}

#endif
