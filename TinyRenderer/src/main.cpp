#include <vector>
#include <iostream>
#include <algorithm>

#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "our_gl.h"

Model *model = NULL;
const int width = 800;
const int height = 800;

Vec3f light_dir{ 1,1,1 };
Vec3f eye{ 1,1,3 };
Vec3f center{ 0,0,0 };
Vec3f up{ 0,1,0 };

struct GouraudShader : public IShader {
	Vec3f varying_intensity; // written by vertex shader, read by fragment shader

	virtual Vec4f vertex(int iface, int nthvert) {
		// get diffuse lighting intensity
		varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert) * light_dir); 
		// read the vertex from .obj file
		Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert));
		//transform it to screen coordinates
		return Viewport * Projection * ModelView * gl_Vertex;
	}

	virtual bool fragment(Vec3f bar, TGAColor &color) {
		// interpolate intensity for the current pixel
		float intensity = varying_intensity * bar;
		// calculate current color
		color = TGAColor(255, 255, 255) * intensity;
		// no, we do not discard this pixel
		return false;
	}
};

struct GShaderMod : public GouraudShader {
	virtual bool fragment(Vec3f bar, TGAColor &color) {
		float intensity = varying_intensity * bar;
		if (intensity > .85) intensity = 1;
		else if (intensity > .60) intensity = .80;
		else if (intensity > .45) intensity = .60;
		else if (intensity > .30) intensity = .45;
		else if (intensity > .15) intensity = .30;
		else intensity = 0;
		color = TGAColor(255, 155, 0) * intensity;
		return false;
	}
};

struct TextureShader : public IShader {
	Vec3f varying_intensity; 	// written by vertex shader, read by fragment shader
	mat<2,3,float> varying_uv;	// same as above 

	virtual Vec4f vertex(int iface, int nthvert) {
		// set 2x3 matrix with vertices idx (0..2) and uv-values (Vec2f)
		varying_uv.set_col(nthvert, model->uv(iface,nthvert));
		// get diffuse lighting intensity
		varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert) * light_dir); 
		// read the vertex from .obj file
		Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert));
		//transform it to screen coordinates
		return Viewport * Projection * ModelView * gl_Vertex;
	}

	virtual bool fragment(Vec3f bar, TGAColor &color) {
		// interpolate intensity for the current pixel
		float intensity = varying_intensity * bar;
		// calculate uv values for the current pixel
		Vec2f uv = varying_uv * bar;
		// calculate current color
		color = model->diffuse(uv) * intensity;
		// no, we do not discard this pixel
		return false;

	}
};

struct Shader : public IShader {
	mat<2,3,float> varying_uv;	// same as above 
	mat<4,4,float> uniform_M;	// Projection*ModelView
	mat<4,4,float> uniform_MIT;	// (Projection*ModelView).invert_transpose()

	virtual Vec4f vertex(int iface, int nthvert) {
		// set 2x3 matrix with vertices idx (0..2) and uv-values (Vec2f)
		varying_uv.set_col(nthvert, model->uv(iface,nthvert));
		// read the vertex from .obj file
		Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert));
		//transform it to screen coordinates
		return Viewport * Projection * ModelView * gl_Vertex;
	}

	virtual bool fragment(Vec3f bar, TGAColor &color) {
		// calculate uv values for the current pixel
		Vec2f uv = varying_uv * bar;
		Vec3f n = proj<3>(uniform_MIT * embed<4>(model->normal(uv))).normalize();
		Vec3f l = proj<3>(uniform_M   * embed<4>(light_dir		  )).normalize();
		// calculate light intensity
		float intensity = std::max(0.f, n * l);
		// calculate current color
		color = model->diffuse(uv) * intensity;
		// no, we do not discard this pixel
		return false;

	}
};

int main(int argc, char** argv) {
	if (2 == argc) {
		model = new Model(argv[1]);
	}
	else {
		model = new Model("obj/african_head.obj");
	}

	lookat(eye, center, up);
	viewport(width/8, height/8, width*3/4, height*3/4);
	projection(-1.f / (eye - center).norm());
	light_dir.normalize();

	TGAImage   image(width, height, TGAImage::RGB);
	TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);

	// GouraudShader shader;
	// GShaderMod shader;
	// TextureShader shader;
	Shader shader;
	shader.uniform_M = Projection * ModelView;
	shader.uniform_MIT = (Projection * ModelView).invert_transpose();

	for (int i = 0; i < model->nfaces(); i++) {
		Vec4f screen_coords[3];
		for (int j = 0; j < 3; j++) {
			screen_coords[j] = shader.vertex(i, j);
		}
		triangle(screen_coords, shader, image, zbuffer);
	}

	image.flip_vertically();
	zbuffer.flip_vertically();
	image.write_tga_file("output.tga");
	zbuffer.write_tga_file("zbuffer.tga");

	delete model;
	return 0;
}