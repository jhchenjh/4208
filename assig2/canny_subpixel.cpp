#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdint.h>
#include <cmath>
#include <cfloat>

#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

int greater_than(float a, float b);

using namespace std;

uint8_t img_header[5];

#define col                     *(uint16_t *)(img_header)
#define row                     *(uint16_t *)(img_header + 2)
#define pixel_size              *(uint8_t *)(img_header + 4)
#define PI 3.14159265

int main(int argc, char **argv)
{
	//Input the raw image, read the raw image to a matrix
	FILE *input;
	uint8_t get_char;
	input = fopen(argv[1],"rb");
	//read the image header
	for (int i = 0; i < 5; i++)
	{
    img_header[i] = fgetc(input);
	}

	//define all image arrays
	int img[row][col]; 			//Original image
	int img_gaussian[row][col]; //Image after gaussian mask
	int img_sobel[row][col];	//Image afte sobel mask
	int img_nms[row][col];  	//Image after non maximum supression

	//read the image in to 2d array
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			img[i][j] = fgetc(input);     //read the image to an array
		}
	}
	printf("Loaded image: %s, no of rows: %d, no of cols: %d\n", argv[1], row, col);

  float sigma = stod(argv[2]);
	int kernal_size;
	if(sigma <0.8){
		kernal_size = 3;
	}
	else if(sigma <1.2){
		kernal_size = 5;
	}
	else if(sigma <1.6){
		kernal_size = 7;
	}
	else if(sigma <2){
		kernal_size = 9;
	}
	else if(sigma <2.4){
		kernal_size = 11;
	}
	else if(sigma <2.8){
		kernal_size = 13;
	}
	else{
		kernal_size = 15;
	}
  printf("Gaussian mask with sigma=%.2f, kernal_size=%d:\n", sigma, kernal_size);

  int gaussianoffset = (kernal_size-1)/2;
	float filter[kernal_size][kernal_size];
	int gaussian_mask[kernal_size][kernal_size];

	for(int x = -gaussianoffset; x <= gaussianoffset; x++)
  {
		for(int y = -gaussianoffset; y <= gaussianoffset; y++)
    {
			float expo =(x*x+y*y)/(2*sigma*sigma);
			filter[x+gaussianoffset][y+gaussianoffset] =  (exp(-expo))/(sigma*sigma*(2*M_PI));
		}
	}

	float temp = filter[0][0];
	int scale = 0;
	for (int i = 0; i < kernal_size; i++)
	{
		for (int j = 0; j < kernal_size; j++)
		{
			gaussian_mask[i][j] = round(filter[i][j] / temp);
			scale = scale + gaussian_mask[i][j];
			printf("%5d   ", gaussian_mask[i][j]);
		}
		printf(" \n");
	}
	printf("The scale is %d \n", scale);

	for(int i = gaussianoffset; i < row-gaussianoffset; i++)
	{
		for(int j = gaussianoffset; j < col-gaussianoffset; j++)
		{
			int sum = 0;
			for(int a = 0; a < kernal_size; a++)
			{
				for(int b = 0; b < kernal_size; b++)
				{
					sum = sum + gaussian_mask[a][b] * img[i-gaussianoffset+a][j-gaussianoffset+b];
				}
			}
			img_gaussian[i][j] = round(sum/scale);
		}
	}
	printf("Image has convolute with gaussian\n");

  int sobelx[3][3] = {{-1, 0, 1},{-2, 0, 2},{-1, 0, 1}};
  int sobely[3][3] = {{-1, -2, -1},{0, 0, 0},{1, 2, 1}};
	int G_x[row][col];  	//Image after sobel mask x
	int G_y[row][col];  	//Image after sobel mask y
  float G_strength[row][col];
  float G_orientation[row][col];

	for (int i = 1; i < row-1; i++)
	{
		for (int j = 1; j < col-1; j++)
		{
			int sum = 0;
			for (int a = 0; a < 3; a++)
			{
				for (int b = 0; b < 3; b++)
				{
					sum = sum + sobelx[a][b]*img_gaussian[i+a-1][j+b-1];
				}
			}
			G_x[i][j] = sum;
		}
	}
	for (int i = 1; i < row-1; i++)
	{
		for (int j = 1; j < col-1; j++)
		{
			int sum = 0;
			for (int a = 0; a < 3; a++)
			{
				for (int b = 0; b < 3; b++)
				{
					sum = sum + sobely[a][b]*img_gaussian[i+a-1][j+b-1];
				}
			}
			G_y[i][j] = sum;
		}
	}
  printf("Image has convolute with gaussian 1\n");

  for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
    {
      G_strength[i][j] = sqrt(pow(G_x[i][j],2)+pow(G_y[i][j],2));
      float angle = atan(float(G_x[i][j])/G_y[i][j])*180/PI;
      G_orientation[i][j] = (angle < 0) ? (angle+180) : angle;
    }
	}

  static float E_x[1000][1000];
  static float E_y[1000][1000];
  int subedge[row][col];

  for(int i = 0; i < row; i++)
	{
		for(int j = 0; j < col; j++)
    {
      E_x[i][j] = -1.0;
      E_y[i][j] = -1.0;
    }
  }

  for(int i=2; i<(row-2); i++)
  {
    for(int j=2; j<(col-2); j++)
    {
      int Dx = 0;
      int Dy = 0;
      float mod = G_strength[i][j];
  		float U = G_strength[i][j-1];
  		float D = G_strength[i][j+1];
  		float L = G_strength[i-1][j];
  		float R = G_strength[i+1][j];
  		float gx = fabs(G_x[i][j]);
  		float gy = fabs(G_y[i][j]);

      if (greater_than(mod, L) && !greater_than(R, mod) && gx >= gy) Dx = 1;
      else if (greater_than(mod, D) && !greater_than(U, mod) && gx <= gy) Dy = 1;

      if (Dx > 0 || Dy > 0)
  		{
  			/* offset value is in [-0.5, 0.5] */
  			float a = G_strength[i-Dx][j-Dy];
  			float b = G_strength[i][j];
  			float c = G_strength[i+Dx][j+Dy];
  			float offset = 0.5 * (a - c) / (a - b - b + c);

  			/* store edge point */
  			E_x[i][j] = i + offset * Dx;
  			E_y[i][j] = j + offset * Dy;
        subedge[i][j] = int(sqrt(pow(E_x[i][j],2)+pow(E_y[i][j],2)));
  		}
      else
      {
        subedge[i][j] = int(G_strength[i][j]);
      }
    }
  }


  // int magnitude_NMS[row][col];
	// for (int i = 1; i < row-1; i++)
	// {
	// 	for (int j = 1; j < col-1; j++)
	// 	{
	// 		if (subedge[i][j] != 0)
	// 		{
	// 			float angle = atan(float(E_y[i][j])/E_x[i][j])*180/PI;
	// 			if (angle < 0)
	// 				angle = angle + 180;
	//
	// 			if (angle < 22.5 || angle >= 157.5)
	// 			{
	// 				if(subedge[i][j]>=subedge[i-1][j] && subedge[i][j]>=subedge[i+1][j])
	// 				{
	// 					magnitude_NMS[i][j] = subedge[i][j];
	// 					magnitude_NMS[i-1][j] = 0;
	// 					magnitude_NMS[i+1][j] = 0;
	// 				}
	// 			}
	//
	// 			if (angle >=22.5 && angle < 67.5)
	// 			{
	// 				if(subedge[i][j]>=subedge[i-1][j-1] && subedge[i][j]>=subedge[i+1][j+1])
	// 				{
	// 					magnitude_NMS[i][j] = subedge[i][j];
	// 					magnitude_NMS[i-1][j-1] = 0;
	// 					magnitude_NMS[i+1][j+1] = 0;
	// 				}
	// 			}
	//
	// 			if (angle >= 67.5 && angle < 112.5)
	// 			{
	// 				if(subedge[i][j]>=subedge[i][j-1] && subedge[i][j]>=subedge[i][j+1])
	// 				{
	// 					magnitude_NMS[i][j] = subedge[i][j];
	// 					magnitude_NMS[i][j-1] = 0;
	// 					magnitude_NMS[i][j+1] = 0;
	// 				}
	// 			}
	//
	// 			if (angle >= 112.5 && angle < 157.5)
	// 			{
	// 				if(subedge[i][j]>=subedge[i-1][j+1] && subedge[i][j]>=subedge[i+1][j-1])
	// 				{
	// 					magnitude_NMS[i][j] = subedge[i][j];
	// 					magnitude_NMS[i-1][j+1] = 0;
	// 					magnitude_NMS[i+1][j-1] = 0;
	// 				}
	// 			}
	//
	// 		}
	// 	}
	// }


  //display img
	FILE* pgmimg;
  pgmimg = fopen(argv[3], "wb");
  fprintf(pgmimg, "P2\n");
  fprintf(pgmimg, "%d %d\n", col, row);
  fprintf(pgmimg, "255\n");
  int count = 0;
  for(int i = 0; i < row; i++)
  {
    for(int j = 0; j < col; j++)
    {
    	int temp = subedge[i][j];
    	if (temp>255)
    		temp = 255;
    	if (temp<0)
    		temp = 0;
        fprintf(pgmimg, "%d ", temp);
    }
    fprintf(pgmimg, "\n");
  }
  fclose(pgmimg);

  return 0;
}


/* compute a > b considering the rounding errors due to the representation of float numbers*/
int greater_than(float a, float b)
{
	if (a <= b) return FALSE;  /* trivial case, return as soon as possible */
	if ((a - b) < 1000 * DBL_EPSILON) return FALSE;
	return TRUE; /* greater_than */
}
