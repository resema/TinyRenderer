#include <vector>
#include <cmath>
#include <algorithm>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
Model *model = NULL;
const int width = 200;
const int height = 200;

//
//  Tutorial 1

//void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color) {
//	bool steep = false;
//	if (std::abs(p0.x - p1.x)<std::abs(p0.y - p1.y)) {
//		std::swap(p0.x, p0.y);
//		std::swap(p1.x, p1.y);
//		steep = true;
//	}
//	if (p0.x>p1.x) {
//		std::swap(p0, p1);
//	}
//
//	for (int x = p0.x; x <= p1.x; x++) {
//		float t = (x - p0.x) / (float)(p1.x - p0.x);
//		int y = p0.y*(1. - t) + p1.y*t;
//		if (steep) {
//			image.set(y, x, color);
//		}
//		else {
//			image.set(x, y, color);
//		}
//	}
//}

// 
//  Tutorial 2

//void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
//	if (t0.y == t1.y && t0.y == t2.y) return;
//	// sort the vertices t0,t1,t2 lower to upper, bubblesort....
//	if (t0.y>t1.y) std::swap(t0, t1);
//	if (t0.y>t2.y) std::swap(t0, t2);
//	if (t1.y>t2.y) std::swap(t1, t2);
//
//	int total_height = t2.y - t0.y;
//	for (int i = 0; i < total_height; i++) {
//		bool second_half = i > t1.y - t0.y || t1.y == t0.y;
//		int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
//		float alpha = (float)i / total_height;
//		float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height;
//		Vec2i A = t0 + (t2 - t0)*alpha;
//		Vec2i B = second_half ? t1 + (t2 - t1)*beta : t0 + (t1 - t0)*beta;
//		if (A.x > B.x) std::swap(A, B);
//		
//		for (int j = A.x; j <= B.x; j++) {
//			image.set(j, t0.y + i, color);
//		}
//	}
//}

//
// Multi threading approach

Vec3f barycenttric(Vec2i *pts, Vec2i P) {
	Vec3f u = cross<float>(Vec3f(pts[2][0] - pts[0][0], pts[1][0] - pts[0][0], pts[0][0] - P[0]), Vec3f(pts[2][1] - pts[0][1], pts[1][1] - pts[0][1], pts[0][1] - P[1]));
	if (std::abs(u[2]) < 1) return Vec3f(-1, 1, 1);
	return Vec3f(1.f - (u.x + u.y) / u.z, u.y/u.z, u.x / u.z);
}

void triangle(Vec2i *pts, TGAImage &image, TGAColor color) {
	Vec2i bboxmin(image.get_width() - 1, image.get_height() - 1);
	Vec2i bboxmax(0, 0);
	Vec2i clamp(image.get_width() - 1, image.get_height() - 1);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			bboxmin[j] = std::max<int>(0, std::min<int>(bboxmin[j], pts[i][j]));
			bboxmax[j] = std::min<int>(clamp[j], std::max<int>(bboxmax[j], pts[i][j]));
		}
	}
	Vec2i P;
	for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
		for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
			Vec3f bc_screen = barycenttric(pts, P);
			if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
			image.set(P.x, P.y, color);
		}
	}
}

int main(int argc, char** argv) {
	if (2 == argc) {
		model = new Model(argv[1]);
	}
	else {
		model = new Model("obj/african_head.obj");
	}
	TGAImage frame(width, height, TGAImage::RGB);
	//Vec2i pts[3] = { Vec2i(10,10),Vec2i(100,30), Vec2i(190,160) };
	//triangle(pts, frame, TGAColor(200, 100, 150, 0));
	for (int i = 0; i < model->nfaces(); i++) {
		std::vector<int> face = model->face(i);
		Vec2i screen_coords[3];
		for (int j = 0; j < 3; j++) {
			Vec3f world_coords = model->vert(face[j]);
			screen_coords[j] = Vec2i((world_coords.x + 1.)*width / 2., (world_coords.y + 1.)*height / 2.);
		}
		triangle(screen_coords, frame, TGAColor(rand() % 255, rand() % 255, rand() % 255, 255));
	}
	
	frame.flip_vertically(); // to place the origin in the bottom left corner of the image
	frame.write_tga_file("framebuffer.tga");

	return 0;
}
//
//  Tutorial 1: Draw face verts

//void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
//	bool steep = false;
//	if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
//		std::swap(x0, y0);
//		std::swap(x1, y1);
//		steep = true;
//	}
//	if (x0 > x1) { // make it left-to-right
//		std::swap(x0, x1);
//		std::swap(y0, y1);
//	}
//
//	for (int x = x0; x <= x1; x++) {
//		float t = (x - x0) / (float)(x1 - x0);
//		int y = y0*(1. - t) + y1*t;
//		if (steep) {
//			image.set(y, x, color);
//		}
//		else {
//			image.set(x, y, color);
//		}
//	}
//}

//int main(int argc, char** argv) {
//	if (2 == argc) {
//		model = new Model(argv[1]);
//	}
//	else {
//		model = new Model("obj/african_head.obj");
//	}
//
//	TGAImage image(width, height, TGAImage::RGB);
//	for (int i = 0; i<model->nfaces(); i++) {
//		std::vector<int> face = model->face(i);
//		for (int j = 0; j<3; j++) {
//			Vec3f v0 = model->vert(face[j]);
//			Vec3f v1 = model->vert(face[(j + 1) % 3]);
//			int x0 = (v0.x + 1.)*width / 2.;
//			int y0 = (v0.y + 1.)*height / 2.;
//			int x1 = (v1.x + 1.)*width / 2.;
//			int y1 = (v1.y + 1.)*height / 2.;
//			line(x0, y0, x1, y1, image, white);
//		}
//	}
//
//	image.flip_vertically(); // to have the origin at the left bottom corner of the image
//	image.write_tga_file("output.tga");
//	delete model;
//	return 0;
//}