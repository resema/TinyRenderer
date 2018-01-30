#include <vector>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <limits>

#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "our_gl.h"

Model *model = NULL;
float *shadowbuffer = NULL;

const int width = 800;
const int height = 800;

//Vec3f light_dir{ 1.f,1.f,1.f };
Vec3f eye{ 1.2f, -0.8f, 3.f };
Vec3f center{ 0,0,0 };
Vec3f up{ 0,1,0 };

TGAImage total(1024, 1024, TGAImage::GRAYSCALE);
TGAImage  occl(1024, 1024, TGAImage::GRAYSCALE);

struct ZShader : public IShader {
	// triangle coordinates (clip coordinates), written by VS, read by FS
	mat<4, 3, float> varying_tri;

	virtual Vec4f vertex(int iface, int nthvert) {
		// read the vertex from .obj file and transform it to screen coordinates
		Vec4f gl_Vertex = Projection * ModelView * embed<4>(model->vert(iface, nthvert));
		varying_tri.set_col(nthvert, gl_Vertex / gl_Vertex[3]);
		// return vertex coordinates
		return gl_Vertex;
	}

	virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor &color) {
		color = TGAColor(0, 0, 0);
		// no, we do not discard this pixel
		return false;
	}
};

float max_elevation_angle(float *zbuffer, Vec2f p, Vec2f dir) {
	float maxangle = 0;
	for (float t = 0.f; t < 1000.f; t += 1.f) {
		Vec2f cur = p + dir*t;
		if (cur.x >= width || cur.y >= height || cur.x < 0 || cur.y < 0)
			return maxangle;

		float distance = (p - cur).norm();
		if (distance < 1.f) continue;
		float elevation = zbuffer[int(cur.x) + int(cur.y) * width] - zbuffer[int(p.x) + int(p.y) * width];
		maxangle = std::max(maxangle, atanf(elevation / distance));
	}

	return maxangle;
}

int main(int argc, char** argv) {
	if (2 > argc) {
		std::cerr << "Usage: " << argv[0] << " obj/model.obj" << std::endl;
        return 1;
	}

	float *zbuffer = new float[width*height];
	for (int i = width*height; --i; ) {
		zbuffer[i] = -std::numeric_limits<float>::max();
	}
	model = new Model(argv[1]);

	TGAImage frame(width, height, TGAImage::RGB);
	lookat(eye, center, up);
	viewport(width/8, height/8, width*3/4, height*3/4);
	projection(-1.f / (eye-center).norm());
	
	ZShader zshader;
	for (int i = 0; i < model->nfaces(); i++) {
		for (int j = 0; j < 3; j++) {
			zshader.vertex(i, j);
		}
		triangle(zshader.varying_tri, zshader, frame, zbuffer);
	}

	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			if (zbuffer[x + y * width] < -1e5) continue;
			float total = 0;
			// 8 ray cast in each direction around the current pixel
			for (float a = 0.f; a < M_PI * 2 - 1e-4; a += M_PI / 4) {
				// solid angle approximation
				total += M_PI / 2 - max_elevation_angle(zbuffer, Vec2f(x, y), Vec2f(cos(a), sin(a)));
			}
			total /= (M_PI / 2) * 8;
			total = pow(total, 100.f);
			frame.set(x, y, TGAColor(total * 255, total * 255, total * 255));
		}
	}
	
	frame.flip_vertically();
	frame.write_tga_file("framebuffer.tga");

	delete model;
	delete [] zbuffer;
	return 0;
}