# Makefile for project raytracing
# RP
#
raytracing: ./source/raytracing.cpp
	g++ -o raytracing ./source/raytracing.cpp -lGL -lGLU -lglut -lGLEW
	mv raytracing ./bin/

clean:
	rm ./bin/*
