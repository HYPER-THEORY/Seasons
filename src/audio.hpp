#ifndef AUDIO_HPP
#define AUDIO_HPP

#include <string>
#include <SDL2/SDL.h>

class audio {
public:
	bool loop = false;
	double volume = 1;
	int position = 0;
	
	static void init() {
		SDL_Init(SDL_INIT_AUDIO);
	}
	
	audio(std::string f, bool l = false) : loop(l) {
		SDL_LoadWAV(f.c_str(), &spec, &buffer, &length);
		spec.userdata = this;
		spec.callback = [](void* ud, uint8_t* s, int l) {
			SDL_memset(s, 0, l);
			audio* data = static_cast<audio*>(ud);
			if (data->position + l < data->length) {
				SDL_MixAudioFormat(s, data->buffer + data->position, data->spec.format, l,
					SDL_MIX_MAXVOLUME * data->volume);
				data->position += l;
			} else {
				SDL_MixAudioFormat(s, data->buffer + data->position, data->spec.format, data->length - data->position,
					SDL_MIX_MAXVOLUME * data->volume);
				data->position = 0;
				if (!data->loop) data->stop();
			}
		};
		device = SDL_OpenAudioDevice(nullptr, 0, &spec, nullptr, 0);
	}
	
	~audio() {
		SDL_CloseAudioDevice(device);
		SDL_FreeWAV(buffer);
	}
	
	void play() {
		SDL_PauseAudioDevice(device, 0);
	}
	
	void stop() {
		SDL_PauseAudioDevice(device, 1);
	}
	
private:
	SDL_AudioSpec spec;
	SDL_AudioDeviceID device;
	uint8_t* buffer = nullptr;
	uint32_t length = 0;
};

#endif
