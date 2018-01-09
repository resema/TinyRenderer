#include "our_gl.h"

void viewport(int x, int y, int w, int h)
{
	
}

void projection(float coeff)
{
	
}

void lookat(Vec3f eye, Vec3f center, Vec3f up)
{
	Vec3f z = (eye - center).normalize();
	Vec3f x = cross(up, z).normalize();
	Vec3f y = cross(z, x).normalize();
	
}


//
// Rasterizer
void triangle(Vec4f *pts, IShader &shader, TGAImage &image, TGAImage &zbuffer)
{

}