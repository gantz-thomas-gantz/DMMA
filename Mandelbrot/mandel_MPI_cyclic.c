/*
 * Sorbonne Université  --  PPAR
 * Computing the Mandelbrot set, sequential version
 */
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>               /* compile with -lm */
#include <sys/time.h>
#include <arpa/inet.h>          /* htonl */

#include "rasterfile.h"

char info[] = "\
Usage:\n\
      ./mandel dimx dimy xmin ymin xmax ymax depth\n\
\n\
      dimx,dimy : dimensions of the image to generate\n\
      xmin,ymin,xmax,ymax : domain to be computed, in the complex plane\n\
      depth : maximal number of iterations\n\
\n\
Some examples of execution:\n\
      ./mandel 3840 3840 0.35 0.355 0.353 0.358 200\n\
      ./mandel 3840 3840 -0.736 -0.184 -0.735 -0.183 500\n\
      ./mandel 3840 3840 -1.48478 0.00006 -1.48440 0.00044 100\n\
      ./mandel 3840 3840 -1.5 -0.1 -1.3 0.1 10000\n\
";

double wallclock_time()
{
	struct timeval tmp_time;
	gettimeofday(&tmp_time, NULL);
	return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6);
}

unsigned char cos_component(int i, double freq)
{
	double iD = i;
	iD = cos(iD / 255.0 * 2 * M_PI * freq);
	iD += 1;
	iD *= 128;
	return iD;
}

/**
 *  Save the data array in rasterfile format
 */
void save_rasterfile(char *name, int largeur, int hauteur, unsigned char *p)
{
	FILE *fd = fopen(name, "w");
	if (fd == NULL) {
		perror("Error while opening output file. ");
		exit(1);
	}

	struct rasterfile file;
	file.ras_magic = htonl(RAS_MAGIC);
	file.ras_width = htonl(largeur);	            /* width of the image, in pixels */
	file.ras_height = htonl(hauteur);	            /* height of the image, in pixels */
	file.ras_depth = htonl(8);	                    /* depth of each pixel (1, 8 or 24 ) */
	file.ras_length = htonl(largeur * hauteur);	    /* size of the image in nb of bytes */
	file.ras_type = htonl(RT_STANDARD);	            /* file type */
	file.ras_maptype = htonl(RMT_EQUAL_RGB);
	file.ras_maplength = htonl(256 * 3);
	fwrite(&file, sizeof(struct rasterfile), 1, fd);

	/* Color palette: red component */
	for (int i = 255; i >= 0; i--) {
		unsigned char o = cos_component(i, 13.0);
		fwrite(&o, sizeof(unsigned char), 1, fd);
	}

	/* Color palette: green component */
	for (int i = 255; i >= 0; i--) {
		unsigned char o = cos_component(i, 5.0);
		fwrite(&o, sizeof(unsigned char), 1, fd);
	}

	/* Color palette: blue component */
	for (int i = 255; i >= 0; i--) {
		unsigned char o = cos_component(i + 10, 7.0);
		fwrite(&o, sizeof(unsigned char), 1, fd);
	}

	fwrite(p, largeur * hauteur, sizeof(unsigned char), fd);
	fclose(fd);
}

/**
 * Given the coordinates of a point $c = a + ib$ in the complex plane,
 * the function returns the corresponding color which estimates the
 * distance from the Mandelbrot set to this point.
 * Consider the complex sequence defined by:
 * \begin{align}
 *     z_0     &= 0 \\
 *     z_{n+1} &= z_n^2 + c
 *   \end{align}
 * The number of iterations that this sequence needs before diverging
 * is the number $n$ for which $|z_n| > 2$.
 * This number is brought back to a value between 0 and 255, thus
 * corresponding to a color in the color palette.
 */
unsigned char xy2color(double a, double b, int depth)
{
	double x = 0;
	double y = 0;
	for (int i = 0; i < depth; i++) {
		/* save the former value of x (which will be erased) */
		double temp = x;
		/* new values for x and y */
		double x2 = x * x;
		double y2 = y * y;
		x = x2 - y2 + a;
		y = 2 * temp * y + b;
		if (x2 + y2 > 4.0)
			return (i % 255);   /* diverges */
	}
	return 255;   /* did not diverge in depth steps */
}

/* 
 * Main: in each point of the grid, apply xy2color
 */
int main(int argc, char **argv)
{
	if (argc == 1) {
		fprintf(stderr, "%s\n", info);
		return 1;
	}

	MPI_Init(NULL, NULL);

	/* Default values for the fractal */
	double xmin = -2;     /* Domain of computation, in the complex plane */
	double ymin = -2;
	double xmax = 2;
	double ymax = 2;
	int w = 3840;         /* dimensions of the image (4K HTDV!) */
	int h = 3840;
	int depth = 10000;     /* depth of iteration */

	/* Retrieving parameters */
	if (argc > 1)
		w = atoi(argv[1]);
	if (argc > 2)
		h = atoi(argv[2]);
	if (argc > 3)
		xmin = atof(argv[3]);
	if (argc > 4)
		ymin = atof(argv[4]);
	if (argc > 5)
		xmax = atof(argv[5]);
	if (argc > 6)
		ymax = atof(argv[6]);
	if (argc > 7)
		depth = atoi(argv[7]);

	/* Computing steps for increments */
	double xinc = (xmax - xmin) / (w - 1);
	double yinc = (ymax - ymin) / (h - 1);

	/* show parameters, for verification
	fprintf(stderr, "Domain: [%g,%g] x [%g,%g]\n", xmin, ymin, xmax, ymax);
	fprintf(stderr, "Increment: %g, %g\n", xinc, yinc);
	fprintf(stderr, "depth: %d\n", depth);
	fprintf(stderr, "Image size: %d x %d\n", w, h);
  */
	/* Allocate memory for the output array */
	unsigned char *image = malloc(w * h);
	if (image == NULL) {
		perror("Error while allocating memory for array: ");
		exit(1);
	}

	/* start timer */
	double start = wallclock_time();

	int my_rank; /* rank of the process */
  int p; /* number of processes */
	  
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	int swp = (w*h)/p;

	/* Process the grid point by point */
	for(int k = 0; k<(h*w); k++){
		if((k%p)==my_rank){
			int j = k % w;
			double x = xmin + j*xinc;
			int i = k / w;
			double y = ymin + i*yinc;
			image[k] = xy2color(x, y, depth);
		}	
	}

  /* stop timer for this process */
  double end = wallclock_time();
  double local_time = end - start;

	unsigned char* reduced_image;


	if(my_rank==0){
		reduced_image = malloc(h*w*sizeof(unsigned char));	
	}

	MPI_Reduce(image, reduced_image, h*w, MPI_UNSIGNED_CHAR, MPI_SUM, 0, MPI_COMM_WORLD);
  
	/* Find the minimum and maximum computation time across all processes */
  double min_time, max_time;
  MPI_Reduce(&local_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	fprintf(stderr, "Rank %d total computing time: %g sec\n", my_rank, end - start);

	if(my_rank==0){

	  fprintf(stderr, "Shortest computation time: %g sec\n", min_time);
    fprintf(stderr, "Longest computation time: %g sec\n", max_time);	
		/* Save the image in the output file "mandel.ras" */
		save_rasterfile("mandel.ras", w, h, image);

	}

	MPI_Finalize();
	return 0;
}

