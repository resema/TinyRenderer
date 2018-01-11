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

Vec3f light_dir{ 1,-1,1 };
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

//Vec3f m2v(Matrix m) {
//	return Vec3f(m[0][0] / m[3][0], m[1][0] / m[3][0], m[2][0] / m[3][0]);
//}
//
//Matrix v2m(Vec3f v) {
//	Matrix m(4, 1);
//	m[0][0] = v.x;
//	m[1][0] = v.y;
//	m[2][0] = v.z;
//	m[3][0] = 1.f;
//	return m;
//}
//
//Matrix modelview(Vec3f eye, Vec3f center, Vec3f up) {
//    Vec3f z = (eye-center).normalize();
//    Vec3f x = (up^z).normalize();
//    Vec3f y = (z^x).normalize();
//    Matrix mv = Matrix::identity(4);
//    for (int i=0; i<3; i++) {
//        mv[0][i] = x[i];
//        mv[1][i] = y[i];
//        mv[2][i] = z[i];
//        mv[i][3] = -center[i];
//    }
//    return mv;
//}

//Matrix viewport(int x, int y, int w, int h) {
//	Matrix m = Matrix::identity(4);
//	m[0][3] = x + w / 2.f;
//	m[1][3] = y + h / 2.f;
//	m[2][3] = depth / 2.f;
//
//	m[0][0] = w / 2.f;
//	m[1][1] = h / 2.f;
//	m[2][2] = depth / 2.f;
//	return m;
//}

//void triangle(Vec3i t0, Vec3i t1, Vec3i t2, Vec2i uv0, Vec2i uv1, Vec2i uv2, float ity0, float ity1, float ity2, TGAImage &image, int *zbuffer) {
//	if (t0.y == t1.y && t0.y == t2.y) return; // i dont care about degenerate triangles
//	if (t0.y > t1.y) { std::swap(t0, t1); std::swap(uv0, uv1); std::swap(ity0, ity1); }
//	if (t0.y > t2.y) { std::swap(t0, t2); std::swap(uv0, uv2); std::swap(ity0, ity2); }
//	if (t1.y > t2.y) { std::swap(t1, t2); std::swap(uv1, uv2); std::swap(ity1, ity2); }
//
//	int total_height = t2.y - t0.y;
//	for (int i = 0; i<total_height; i++) {
//		bool second_half = i>t1.y - t0.y || t1.y == t0.y;
//		int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
//		float alpha = (float)i / total_height;
//		float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height; // be careful: with above conditions no division by zero here
//		Vec3i A = t0 + Vec3f(t2 - t0)*alpha;
//		Vec3i B = second_half ? t1 + Vec3f(t2 - t1)*beta : t0 + Vec3f(t1 - t0)*beta;
//		Vec2i uvA = uv0 + (uv2 - uv0)*alpha;
//		Vec2i uvB = second_half ? uv1 + (uv2 - uv1)*beta : uv0 + (uv1 - uv0)*beta;
//		float ityA = ity0 + (ity2 - ity0)*alpha;
//		float ityB = second_half ? ity1 + (ity2 - ity1)*beta : ity0 + (ity1 - ity0)*beta;
//		if (A.x > B.x) { std::swap(A, B); std::swap(uvA, uvB); std::swap(ityA, ityB); }
//		for (int j = A.x; j <= B.x; j++) {
//			float phi = B.x == A.x ? 1. : (float)(j - A.x) / (float)(B.x - A.x);
//			Vec3i   P = Vec3f(A) + Vec3f(B - A)*phi;
//			Vec2i uvP = uvA + (uvB - uvA)*phi;
//			float ityP = ityA + (ityB - ityA)*phi;
//			int idx = P.x + P.y*width;
//			if (P.x >= width || P.y >= height || P.x < 0 || P.y < 0) continue;
//			if (zbuffer[idx]<P.z) {
//				zbuffer[idx] = P.z;
//				TGAColor color = model->diffuse(uvP);
//				image.set(P.x, P.y, TGAColor(color.bgra[2], color.bgra[1], color.bgra[0])*ityP);
//			}
//		}
//	}
//}

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

	GouraudShader shader;
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