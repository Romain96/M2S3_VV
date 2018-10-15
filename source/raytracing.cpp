// main.cpp : Test different raytracing algorithms.
//

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../include/scene.h"

const int RX=800;
const int RY=800;

scene<float> sc(RX,RY);


// CALLBACKS for GLUT window
void display() { sc.display(); }
void clic(int button, int state, int x, int y) { sc.clic(button,state,x,y); }
void bouge(int x, int y) { sc.bouge(x,y); }
void key(unsigned char c, int x, int y) { sc.key(c,x,y); }

//  PROGRAMM LOADS A SCENE and WAITS for EVENTS
int main(int argc, char **argv)
{
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowPosition(0,0) ;
  glutInitWindowSize(RX, RY);
  int win = glutCreateWindow("ray tracing");
  GLenum err = glewInit();
  if (GLEW_OK != err)
	{
	/* Problem: glewInit failed, something is seriously wrong. */
	fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
	}

	// init callbacks for interaction
  glutDisplayFunc(display);
  glutMouseFunc(clic);
  glutMotionFunc(bouge);
  glutKeyboardFunc(key);

  // import the scene: .obj format
 sc.importOBJ("scene");

  // classical ray tracing
  sc.pushLight(vec4<float>(0.0f, 10.0f, 9.5f, 1.0f), luminance(1.0f));
  

  // load shaders and objects in OpenGL
  sc.initGL();

  // enter callback - wait for events

  printf("\nUSAGE:\n\nleft mouse button for rotation\n");
  printf("right mouse button for zoom\n");
  printf("middle mouse button for look at\n\n");
  printf("i key for ray tracing a picture\n");
  printf("r key for tracking rays\n");
  printf("w key for switching wire frame / filled display mode\n");
  glutMainLoop();
  return 0;            
}


