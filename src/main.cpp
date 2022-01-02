#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <SDL2/SDL.h>
#include "audio.hpp"
#include "flat.hpp"
#include "gl.hpp"
#include "leaf.hpp"
#include "obj.hpp"
#include "physics.hpp"
#include "window.hpp"
#include "game.hpp"

// core parts
solid character;
flatcamera maincamera;
flatshader natureshader;
std::vector<light*> lights;
std::vector<instance*> instances;
std::vector<solid> solids;
std::vector<std::vector<size_t>> groups;

// audio part
std::vector<audio*> audios;

// particle part
std::vector<leaf> leaves;

// model part
objmodel* block1 = new objmodel();
objmodel* block2 = new objmodel();
objmodel* block3 = new objmodel();
objmodel* tree1 = new objmodel();
objmodel* tree2 = new objmodel();
objmodel* tree3 = new objmodel();
objmodel* tree4 = new objmodel();
objmodel* square1 = new objmodel();
objmodel* square2 = new objmodel();
objmodel* square3 = new objmodel();

// camera part
double axisy = 0;
double axisz = 0;
double charactergradient = 0;
double targetgradient = 0;
vec3 characterdelta = {0, .2, 0};
vec3 targetdelta = {0, .2, 0};
vec3 position = {0, 0, 0};
vec3 direction = {0, 0, 1};
vec3 up = {0, 1, 0};

// camera parameters
uint32_t backgroundcolor = 0xa9d9cb;
double mousesensitivity = .001;
double gradientacceleration = .25;
double deltaacceleration = .2;
vec3 characterbox = {.2, .5, .2};

// kinematics part
bool grounded = false;
bool lastgrounded = false;
int leaveground = 0;
int jumpaction = -1;
int dashaction = -1;
int dashbuffer = 30;
vec3 velocity = {0, 0, 0};

// kinematics parameters
double speed = .075;
double friction = .05;
double collisionloss = .75;
vec3 gravity = {0, -0.01, 0};
vec3 jumpforce = {0, .15, 0};
vec3 limitvelocity = {1, 1, 1};

void solidify(const std::vector<instance*>& is, const std::vector<size_t>& g, solid& s) {
	vec3 v1 = vec3{1, 1, 1} * std::numeric_limits<double>::max();
	vec3 v2 = vec3{1, 1, 1} * -std::numeric_limits<double>::max();
	mat<4, 4> translation;
	mat<4, 4> rotation;
	mat<4, 4> scaling;
	for (auto& i : g) {
		camera::modeltransform(*is[i], translation, rotation, scaling);
		mat<4, 4> affine = translation * rotation * scaling;
		for (auto& v : is[i]->data->vertices) {
			mat<4, 1> matrix = {v.x, v.y, v.z, 1};
			matrix = affine * matrix;
			v1.x = fmin(v1.x, matrix[0][0]);
			v1.y = fmin(v1.y, matrix[1][0]);
			v1.z = fmin(v1.z, matrix[2][0]);
			v2.x = fmax(v2.x, matrix[0][0]);
			v2.y = fmax(v2.y, matrix[1][0]);
			v2.z = fmax(v2.z, matrix[2][0]);
		}
	}
	s = solid(v1, v2.x - v1.x, v2.y - v1.y, v2.z - v1.z);
}

void platform(model* m1, model* m2, model* m3, const vec3& p, double w, double h, double d, double t1, double t2,
	std::vector<instance*>& is, std::vector<std::vector<size_t>>& gs) {
	size_t index = is.size();
	is.push_back(new instance(m1, {p.x, p.y - t1, p.z}, {}, {w, t1, d}));
	is.push_back(new instance(m2, {p.x, p.y - t1 - t2, p.z}, {}, {w, t2, d}));
	if (h - t1 - t2 > eps) {
		is.push_back(new instance(m3, {p.x, p.y - h, p.z}, {}, {w, h - t1 - t2, d}));
		gs.push_back({index, index + 1, index + 2});
	} else {
		gs.push_back({index, index + 1});
	}
}

void tree(model* m, const vec3& p, const vec3& r, const vec3& s, std::vector<instance*>& is,
	std::vector<std::vector<size_t>>& gs) {
	gs.push_back({is.size()});
	is.push_back(new instance(m, p, r, s));
}

void scene(const std::string& s, std::vector<instance*>& is, std::vector<std::vector<size_t>>& gs) {
	// load models
	readobj(s + "block_light.obj", *block1);
	readobj(s + "block.obj", *block2);
	readobj(s + "block_dark.obj", *block3);
	readobj(s + "tree_blocks.obj", *tree1);
	readobj(s + "tree_fat.obj", *tree2);
	readobj(s + "tree_small.obj", *tree3);
	readobj(s + "tree_thin.obj", *tree4);
	readobj(s + "square_light.obj", *square1);
	readobj(s + "square.obj", *square2);
	readobj(s + "square_dark.obj", *square3);
	// platform 1
	platform(block1, block2, block3, {0, 10, 5}, 4, 200, 14, .2, 1, is, gs);
	tree(tree1, {-1, 10, 0}, {}, {1, 1, 1}, is, gs);
	tree(tree2, {1.2, 10, 3}, {}, {1, 1, 1}, is, gs);
	tree(tree3, {-0.6, 10, 5}, {}, {1, 1, 1}, is, gs);
	tree(tree4, {.8, 10, 7.5}, {}, {1, 1, 1}, is, gs);
	tree(tree1, {-1.2, 10, 10}, {}, {1, 1, 1}, is, gs);
	tree(tree2, {1, 10, 10.5}, {}, {1, 1, 1}, is, gs);
	tree(tree3, {1.5, 10, -1.2}, {}, {1, 1, 1}, is, gs);
	is.push_back(new instance(square2, {1, 11.5, 3}, {M_PI_2, M_PI * .025, 0}, {200, 1, .8}));
	is.push_back(new instance(square3, {0, 11.5, 6}, {M_PI_2, -M_PI * .03, 0}, {200, 1, .8}));
	is.push_back(new instance(square3, {1, 12, 6}, {M_PI_2, M_PI * .25, 0}, {200, 1, .4}));
	is.push_back(new instance(square1, {-1, 11.5, 9}, {M_PI_2, M_PI * .02, 0}, {200, 1, .4}));
	is.push_back(new instance(square1, {-1, 12, 9}, {M_PI_2, M_PI * .02, 0}, {200, 1, .4}));
	// platform 2
	platform(block1, block2, block3, {0, 10, 24}, 4, 200, 4, .2, 1, is, gs);
	tree(tree1, {-1, 10, 25}, {}, {1, 1, 1}, is, gs);
	tree(tree2, {-1, 10, 22.5}, {}, {1, 1, 1}, is, gs);
	tree(tree3, {1.5, 10, 23.5}, {}, {1, 1, 1}, is, gs);
	// path between platform 1 and 2
	platform(block1, block2, block3, {0, 10, 14}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {.4, 10, 16}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {-0.4, 10, 18}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {0, 10, 20}, 1, .4, 1, .2, .2, is, gs);
	// platform 3
	platform(block1, block2, block3, {-10, 10, 34}, 4, 200, 4, .2, 1, is, gs);
	tree(tree2, {-11, 10, 35}, {}, {1, 1, 1}, is, gs);
	tree(tree3, {-10, 10, 33}, {}, {1, 1, 1}, is, gs);
	tree(tree4, {-9, 10, 34.8}, {}, {1, 1, 1}, is, gs);
	// path between platform 2 and 3
	platform(block1, block2, block3, {0, 9.6, 28}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {0, 10, 30}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {0, 9.6, 32}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {0, 10, 34}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {-2, 9.6, 34}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {-4, 10, 34}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {-6, 9.6, 34}, 1, .4, 1, .2, .2, is, gs);
	// platform 4
	platform(block1, block2, block3, {-16, 15.4, 34}, 4, 200, 4, .2, 1, is, gs);
	tree(tree1, {-15, 15.4, 33}, {}, {1, 1, 1}, is, gs);
	tree(tree2, {-15, 15.4, 35}, {}, {1, 1, 1}, is, gs);
	tree(tree3, {-17, 15.4, 33}, {}, {1, 1, 1}, is, gs);
	tree(tree4, {-17, 15.4, 35}, {}, {1, 1, 1}, is, gs);
	// path between platform 3 and 4
	platform(block1, block2, block3, {-13.5, 10, 34}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {-13.5, 10.6, 36.5}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {-16, 11.2, 36.5}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {-18.5, 11.8, 36.5}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {-18.5, 12.4, 34}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {-18.5, 13, 31.5}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {-16, 13.6, 31.5}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {-13.5, 14.2, 31.5}, 1, .4, 1, .2, .2, is, gs);
	platform(block1, block2, block3, {-13.5, 14.8, 34}, 1, .4, 1, .2, .2, is, gs);
	// secret platform 1
	platform(block1, block2, block3, {0, 10, -16}, 3, 1.2, 3, .2, 1, is, gs);
	tree(tree1, {0, 10, -16}, {}, {1, 1, 1}, is, gs);
	// background 1
	platform(block1, block2, block3, {-30, 10, -10}, 2, 1.2, 2, .2, 1, is, gs);
	tree(tree1, {-30, 10, -10}, {}, {1, 1, 1}, is, gs);
	// background 2
	platform(block1, block2, block3, {20, 8, 15}, 2, 1.2, 2, .2, 1, is, gs);
	tree(tree2, {20, 8, 15}, {}, {1, 1, 1}, is, gs);
	// background 3
	platform(block1, block2, block3, {-20, 10, 20}, 2, 1.2, 2, .2, 1, is, gs);
	tree(tree3, {-20, 10, 20}, {}, {1, 1, 1}, is, gs);
	// background 4
	platform(block1, block2, block3, {20, 10, -25}, 2, 1.2, 2, .2, 1, is, gs);
	tree(tree4, {20, 10, -25}, {}, {1, 1, 1}, is, gs);
	// background 5
	platform(block1, block2, block3, {-35, 12, 10}, 2, 1.2, 2, .2, 1, is, gs);
	tree(tree1, {-35, 12, 10}, {}, {1, 1, 1}, is, gs);
	// background 6
	platform(block1, block2, block3, {20, 8, -5}, 2, 1.2, 2, .2, 1, is, gs);
	tree(tree2, {20, 8, -5}, {}, {1, 1, 1}, is, gs);
	// background 7
	platform(block1, block2, block3, {-20, 12, -20}, 2, 1.2, 2, .2, 1, is, gs);
	tree(tree3, {-20, 12, -20}, {}, {1, 1, 1}, is, gs);
	// background 8
	platform(block1, block2, block3, {15, 10, 25}, 2, 1.2, 2, .2, 1, is, gs);
	tree(tree4, {15, 10, 25}, {}, {1, 1, 1}, is, gs);
}

void music(const std::string& s, std::vector<audio*>& as) {
	audio::init();
	as.push_back(new audio(s + "smile_sower.wav", true));
	as.push_back(new audio(s + "glass_004.wav"));
	as[0]->play();
}

double decelerate(double v, double f) {
	if (v > 0) return fmax(v - f, 0);
	if (v < 0) return fmin(v + f, 0);
	return v;
}

void restart() {
	character.position = characterbox * -0.5 + vec3{0, 10.25, 0};
//	character.position = characterbox * -0.5 + vec3{-10, 10.25, 34};
	maincamera.loadtime = 0;
	axisy = 0;
	axisz = 0;
}

void load() {
	// initialize music
//	music("/Users/hypertheory/GL/audios/", audios);
	// initialize scene
	scene("/Users/hypertheory/GL/models/", instances, groups);
	solids = std::vector<solid>(instances.size());
	for (size_t i = groups.size(); i --> 0;) solidify(instances, groups[i], solids[i]);
	// initialize camera
	character = solid({}, characterbox.x, characterbox.y, characterbox.z);
	maincamera = flatcamera(width, height, M_PI / 180 * 75, .05, 1000);
	maincamera.type = "none";
	std::copy_n(std::vector<double>{
		0, .5, .8, .86, .89, .91, .93, .94, .95, .96, .97, .98, .985, .99, .995, .1,
	}.data(), 16, maincamera.distances);
	maincamera.cursor[0] = {
		1, 1, 0, 0, 0, 0, 0, 0, 1, 1,
		1, 1, 0, 1, 1, 1, 1, 0, 1, 1,
		0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
		0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
	};
	maincamera.cursor[1] = {
		1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
		0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
		1, 0, 0, 1, 1, 1, 1, 0, 0, 1,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	};
	restart();
	// initalize particle
	for (int i = 0; i < 500; ++i) {
		leaf newleaf(new instance(square1), 10, .005, {.05, .05, .05}, {0, -0.005, 0});
		newleaf.reset(character.position);
		leaves.push_back(newleaf);
		instances.push_back(newleaf.graphic);
	}
}

void update(double dt) {
	// particles
	for (auto& l : leaves) l.update(character.position);
	// mouse event
	axisy -= (window::mousex - width * scale / 2) * mousesensitivity;
	axisz -= (window::mousey - height * scale / 2) * mousesensitivity;
	if (axisz > M_PI_2) axisz = M_PI_2;
	if (axisz < -M_PI_2) axisz = -M_PI_2;
	direction = {sin(axisy) * cos(axisz), sin(axisz), cos(axisy) * cos(axisz)};
	up = {-sin(axisy) * sin(axisz), cos(axisz), -cos(axisy) * sin(axisz)};
	// keyboard event
	int movex = 0;
	int movez = 0;
	vec3 force = gravity;
	if (window::keydown[SDLK_w]) ++movez;
	if (window::keydown[SDLK_s]) --movez;
	if (window::keydown[SDLK_d]) ++movex;
	if (window::keydown[SDLK_a]) --movex;
	if (window::keypressed[SDLK_SPACE]) jumpaction = 3;
	if (window::keypressed[1] || window::keypressed[SDLK_TAB]) dashaction = 3;
	if (window::keypressed[SDLK_r]) restart();
	// correct camera gradient
	targetgradient = movex * M_PI * .02;
	charactergradient += (targetgradient - charactergradient) * gradientacceleration;
	up = up.rotate(direction, charactergradient);
	// correct camera position
	characterdelta.x += (targetdelta.x - characterdelta.x) * deltaacceleration;
	characterdelta.y += (targetdelta.y - characterdelta.y) * deltaacceleration;
	characterdelta.z += (targetdelta.z - characterdelta.z) * deltaacceleration;
	// jump action
	++leaveground;
	if (grounded) leaveground = 0;
	if (jumpaction >= 0) {
		--jumpaction;
		if (grounded || (leaveground <= 3 && velocity.y <= 0)) {
			grounded = false;
			force += jumpforce;
		}
	}
	// dash action
	if (dashaction >= 0) {
		--dashaction;
		if (dashbuffer == 45) {
			velocity *= .5;
			force += direction * .2;
			dashbuffer = 0;
			maincamera.dashtime = 0;
		}
	}
	if (dashbuffer < 45) {
		maincamera.cursortype = 1;
		++dashbuffer;
	} else {
		maincamera.cursortype = 0;
	}
	// apply force to velocity
	vec3 lastposition = character.position;
	if (fabs(velocity.x + force.x) < limitvelocity.x) velocity.x += force.x;
	if (fabs(velocity.y + force.y) < limitvelocity.y) velocity.y += force.y;
	if (fabs(velocity.z + force.z) < limitvelocity.z) velocity.z += force.z;
	if (grounded) {
		velocity.x = decelerate(velocity.x, friction);
		velocity.z = decelerate(velocity.z, friction);
	};
	// character movement
	vec3 move(movex, 0, movez);
	if (movex != 0 || movez != 0) {
		move = move.normalize() * speed;
		move = {move.z * sin(axisy) - move.x * cos(axisy), 0, move.z * cos(axisy) + move.x * sin(axisy)};
	}
	character.move(move + velocity);
	// collision detection
	grounded = false;
	if (character.position.y < lastposition.y + velocity.y - eps) {
		velocity.y = -velocity.y * (1 - collisionloss);
	} else if (character.position.y > lastposition.y + velocity.y + eps) {
		if (!lastgrounded) characterdelta.y += velocity.y;
		grounded = true;
		velocity.y = 0;
	}
	lastgrounded = grounded;
	if (character.position.x < lastposition.x + velocity.x + move.x - eps ||
		character.position.x > lastposition.x + velocity.x + move.x + eps)
		velocity.x = -velocity.x * (1 - collisionloss);
	if (character.position.z < lastposition.z + velocity.z + move.z - eps ||
		character.position.z > lastposition.z + velocity.z + move.z + eps)
		velocity.z = -velocity.z * (1 - collisionloss);
	// set camera extrinsic parameters
	position = character.position + characterbox * .5 + characterdelta;
	maincamera.lookat(position, direction, up);
	// output render result
	for (auto i = instances.begin(); i != instances.end();) *i != nullptr ? ++i : i = instances.erase(i);
	sort(instances.begin(), instances.end(), [](instance* a, instance* b) {
		return position.distance(a->position) < position.distance(b->position);
	});
	fill(framebuffer.begin(), framebuffer.end(), backgroundcolor);
	maincamera.draw(instances, natureshader, framebuffer);
}

int main(int argc, char** argv) {
	width = 960;
	height = 540;
	scale = 1;
	title = "Seasons";
	fps = 30;
	hidecursor = true;
	firstperson = true;
	run();
}
