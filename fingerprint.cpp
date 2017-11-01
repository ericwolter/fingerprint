#include <stdio.h>
#include <malloc.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "stb_image_resize.h"

#define IMAGE_SIZE 256
#define M_LOG2E 1.44269504088896340736 // log2(e)
#define img(r,c) (img[(r)*width + (c)])
#define section(r,c) (section[(r)*stride + (c)])

void printImage(const float* img, const int width, const int height) {
  for (int r = 0; r < height; r++) {
    for (int c = 0; c < width; c++) {
      printf("%.5f ", img(r,c));
    }
    printf("\n");
  }
  printf("\n");
}

float* grayscale(const unsigned char* colorImage, const int width, const int height, const int comp) {
  // printf("Converting to grayscale...\n");
  float *grayImage = (float*)malloc(sizeof(float) * width * height);

  unsigned char red, green, blue = 0;
  for (int r = 0; r < height; r++) {
    for (int c = 0; c < width; c++) {
      red = colorImage[r * width + c * comp + 1];
      green = colorImage[r * width + c * comp + 2];
      blue = colorImage[r * width + c * comp + 3];

      // luminosity method
      grayImage[r * width + c] = 0.21f * red + 0.72f * green + 0.07f * blue;
    }
  }

  return grayImage;
}

float* reduceSize(const float *image, const int width, const int height) {
  // printf("Reducing size to %ix%i...\n", IMAGE_SIZE, IMAGE_SIZE);
  float *scale = (float*)malloc(sizeof(float) * width * height);
  stbir_resize_float(image, width, height, 0, scale, IMAGE_SIZE, IMAGE_SIZE, 0 , 1);
  return scale;
}

void computeHaarWavelet(float *img, const int width, const int stride) {

  // printf("Row by row...\n");
  float *section = (float*)malloc(sizeof(float) * stride * stride);
  for (int r = 0; r < stride; r++) {
    for (int c = 0; c < stride; c++) {
      section(r,c) = img(r,c);
    }
  }

  float average = 0.0;
  float difference = 0.0;
  for (int r = 0; r < stride; r++) {
    for (int c = 0; c < stride; c+=2) {
      const int atr = r;
      const int atc = c / 2;
      const int dtr = r;
      const int dtc = stride / 2 + c / 2;

      // printf("compute average at %i:%i=%f and %i:%i=%f and update at %i:%i\n",r,c,img(r, c),r,c+1,img(r, c + 1),atr,atc);
      // printf("compute difference at %i:%i=%f and %i:%i=%f and update at %i:%i\n",r,c,img(r, c),r,c+1,img(r, c + 1), dtr, dtc);
      average = (section(r, c) + section(r, c + 1));
      difference = (section(r, c) - section(r, c + 1));

      img(atr, atc) = average;
      img(dtr, dtc) = difference;
    }
    // printImage(img, width, height);
  }

  // printf("Column by column...\n");
  for (int r = 0; r < stride; r++) {
    for (int c = 0; c < stride; c++) {
      section(r,c) = img(r,c);
    }
  }

  for (int c = 0; c < stride; c++) {
    for (int r = 0; r < stride; r+=2) {
      const int atr = r / 2;
      const int atc = c;
      const int dtr = stride / 2 + r / 2;
      const int dtc = c;

      // printf("compute average at %i:%i=%f and %i:%i=%f and update at %i:%i\n",r,c,img(r, c),r,c+1,img(r + 1, c),atr,atc);
      // printf("compute difference at %i:%i=%f and %i:%i=%f and update at %i:%i\n",r,c,img(r, c),r,c+1,img(r + 1, c), dtr, dtc);
      average = (section(r, c) + section(r + 1, c)) / 2.0f;
      difference = (section(r, c) - section(r + 1, c)) / 2.0f;

      img(atr, atc) = average;
      img(dtr, dtc) = difference;
    }
    // printImage(img, width, height);
  }

  free(section);
}

float* createHaarImage(float *img, const int width, const int height, const size_t size) {
  // printf("Creating haar image...\n", size, size);
  float *haar = (float*)malloc(sizeof(float) * width * height);
  memcpy(haar, img, sizeof(float) * width * height);

  // TODO: Check if size is power of 2
  int max_level = (int)(log(size) * M_LOG2E);
  // printf("Max level: %i\n", max_level);
  for (int l = 0; l < max_level; l++) {
    int stride = (int)(size / pow(2, l));
    // printf("Computing level %i with stride %i\n", l, stride);
    computeHaarWavelet(haar, width, stride);
  }

  return haar;
}

float* getFingerprint(const char* filepath, int *err) {
  // printf("Computing fingerprint for %s\n", filepath);
  *err = 0;

  int width, height;
  const int desired_channels = 3;
  unsigned char *colorImage = stbi_load(filepath, &width, &height, NULL, 3);
  if (!colorImage) {
    fprintf(stderr, "image loading error: %s", stbi_failure_reason());
    *err = 1;
    return NULL;
  }

  float *grayImage = grayscale(colorImage, width, height, desired_channels);
  stbi_image_free(colorImage);

  float *scaledImage = reduceSize(grayImage, width, height);
  free(grayImage);

  float *haar = createHaarImage(scaledImage, IMAGE_SIZE, IMAGE_SIZE, IMAGE_SIZE);
  free(scaledImage);

  float *fingerprint = (float*)malloc(sizeof(float) * 4 * 4);
  for (size_t r = 0; r < 4; r++) {
    for (size_t c = 0; c < 4; c++) {
      fingerprint[(r)*4 + (c)] = haar[(r)*IMAGE_SIZE + (c)];
    }
  }
  free(haar);

  // normalize fingerprint
  float old_min = 255;
  float old_max = -255;
  float f = 0;
  for (size_t i = 0; i < 16; i++) {
    f = fingerprint[i];
    if (f > old_max) {
      old_max = f;
    }
    if (f < old_min) {
      old_min = f;
    }
  }
  float old_range = old_max - old_min;
  for (size_t i = 0; i < 16; i++) {
    fingerprint[i] = (fingerprint[i] - old_min) / old_range;
  }
  // printImage(fingerprint, 16, 1);

  return fingerprint;
}

int main(int argc, char *argv[]) {
  if(argc < 3) {
    printf("Usage: %s <image1.jpg> <image2.jpg>\n", argv[0]);
    return 2;
  }

  // FILE *f;
  // f = fopen(argv[1], "rb");
  // f = fopen(argv[2], "rb");

  int err;
  float *fingerprint1 = getFingerprint(argv[1], &err);
  if(err) return err;
  float *fingerprint2 = getFingerprint(argv[2], &err);
  if(err) return err;

  float diff_l1 = 0;
  for (size_t i = 0; i < 16; i++) {
    diff_l1 += fabsf(fingerprint1[i] - fingerprint2[i]);
  }
  free(fingerprint1);
  free(fingerprint2);
  diff_l1 /= 16;
  diff_l1 = 1 - diff_l1;
  diff_l1 *= 100;

  printf("%.5f",diff_l1);

  return 0;
}
