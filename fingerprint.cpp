#include <malloc.h>
#include <math.h>
#include <stdio.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "stb_image_resize.h"

#define M_LOG2E 1.44269504088896340736 // log2(e)
#define img(r, c) (img[(r)*width + (c)])
#define section(r, c) (section[(r)*stride + (c)])

void printImage(const double *img, const int width, const int height) {
  for (int r = 0; r < height; r++) {
    for (int c = 0; c < width; c++) {
      printf("%.5f ", img(r, c));
    }
    printf("\n");
  }
  printf("\n");
}

void printHash(const unsigned char *hash, const int size) {
  for (int i = 0; i < size; i++) {
    printf("%d", hash[i]);
  }
  printf("\n");
}

double luminosity(const int red, const int green, const int blue) {
  // based on ITU-R 601-2
  return 299.0 / 1000 * red + 587.0 / 1000 * green + 114.0 / 1000 * blue;
}

double *grayscale(const unsigned char *colorImage, const int width,
                  const int height, const int comp) {
  // printf("Converting to grayscale...\n");
  double *grayImage = (double *)malloc(sizeof(double) * width * height);

  unsigned char red, green, blue;
  for (int r = 0; r < height; r++) {
    for (int c = 0; c < width; c++) {
      red = colorImage[r * width + c * comp + 0];
      green = colorImage[r * width + c * comp + 1];
      blue = colorImage[r * width + c * comp + 2];

      grayImage[r * width + c] = floor(luminosity(red, green, blue)) / 255.0;
    }
  }

  return grayImage;
}

unsigned char *reduceSize(const unsigned char *image, const int width,
                          const int height, const int channels,
                          const int scale) {
  unsigned char *scaled =
      (unsigned char *)malloc(sizeof(unsigned char) * scale * scale * channels);
  int err = stbir_resize_uint8(image, width, height, 0, scaled, scale, scale, 0,
                               channels);
  if (err != 1) {
    fprintf(stderr, "image scaling error: %d", err);
  }
  return scaled;
}

void computeHaarWavelet(double *img, const int width, const int stride) {
  double *section = (double *)malloc(sizeof(double) * stride * stride);
  for (int r = 0; r < stride; r++) {
    for (int c = 0; c < stride; c++) {
      section(r, c) = img(r, c);
    }
  }

  double average = 0.0;
  double difference = 0.0;
  for (int r = 0; r < stride; r++) {
    for (int c = 0; c < stride; c += 2) {
      const int atr = r;
      const int atc = c / 2;
      const int dtr = r;
      const int dtc = stride / 2 + c / 2;

      average = (section(r, c) + section(r, c + 1));
      difference = (section(r, c) - section(r, c + 1));

      img(atr, atc) = average;
      img(dtr, dtc) = difference;
    }
  }

  for (int r = 0; r < stride; r++) {
    for (int c = 0; c < stride; c++) {
      section(r, c) = img(r, c);
    }
  }

  for (int c = 0; c < stride; c++) {
    for (int r = 0; r < stride; r += 2) {
      const int atr = r / 2;
      const int atc = c;
      const int dtr = stride / 2 + r / 2;
      const int dtc = c;

      average = (section(r, c) + section(r + 1, c)) / 2.0f;
      difference = (section(r, c) - section(r + 1, c)) / 2.0f;

      img(atr, atc) = average;
      img(dtr, dtc) = difference;
    }
  }

  free(section);
}

double *createHaarImage(double *img, const int width, const int height,
                        const int level) {
  double *haar = (double *)malloc(sizeof(double) * width * height);
  memcpy(haar, img, sizeof(double) * width * height);

  // TODO: Check if size is power of 2
  int size = width;
  // # 1 -> 256 -> 128x128
  // # 2 -> 128 ->  64x64
  // # 3 ->  64 ->  32x32
  // # 4 ->  32 ->  16x16
  // # 5 ->  16 ->   8x8
  // # 6 ->   8 ->   4x4
  // # 7 ->   4 ->   2x2
  // # 8 ->   2 ->   1x1
  for (int l = 0; l < level; l++) {
    int stride = (int)(size / pow(2, l) + 1e-9);
    computeHaarWavelet(haar, width, stride);
  }

  return haar;
}

unsigned char *loadRGBImage(int *err, const char *filepath, int *width,
                            int *height, int *channels) {
  *channels = 3;
  unsigned char *image = stbi_load(filepath, width, height, NULL, *channels);
  if (!image) {
    fprintf(stderr, "image loading error: %s", stbi_failure_reason());
    *err = 1;
    return NULL;
  }
  return image;
}

int cmpfunc(const void *a, const void *b) {
  double va = *(const double *)a;
  double vb = *(const double *)b;
  return (-1) * ((va > vb) - (va < vb));
}

double *truncateCoefficients(const double *haar, const int scale,
                             const int nCoeffs) {

  const int count = scale * scale;

  double *sorted = (double *)malloc(sizeof(double) * count);
  for (int i = 0; i < count; i++) {
    sorted[i] = fabs(haar[i]);
  }
  qsort(sorted, count, sizeof(double), cmpfunc);

  double *truncated = (double *)malloc(sizeof(double) * count);
  double coefficent;
  for (int i = 0; i < count; i++) {
    coefficent = haar[i];
    if (fabs(coefficent) > sorted[nCoeffs]) {
      truncated[i] = coefficent;
    } else {
      truncated[i] = 0;
    }
  }

  free(sorted);

  return truncated;
}

unsigned char *quantifyMedian(double *fingerprint, const int hashSize) {
  double *sorted = (double *)malloc(sizeof(double) * hashSize);
  for (int i = 0; i < hashSize; i++) {
    sorted[i] = fingerprint[i];
  }
  qsort(sorted, hashSize, sizeof(double), cmpfunc);

  unsigned char *quantified =
      (unsigned char *)malloc(sizeof(unsigned char) * hashSize);

  double median;
  if (hashSize % 2 == 0) {
    median = ((sorted[hashSize / 2] + sorted[hashSize / 2 - 1]) / 2.0);
  } else {
    median = sorted[hashSize / 2];
  }

  for (int i = 0; i < hashSize; i++) {
    if (fingerprint[i] > median) {
      quantified[i] = 1;
    } else {
      quantified[i] = 0;
    }
  }

  return quantified;
}

unsigned char *getFingerprint(int *err, const char *filepath, int hashSize) {

  // 0. Image
  int width, height, channels;
  unsigned char *colorImage =
      loadRGBImage(err, filepath, &width, &height, &channels);
  if (*err) {
    return NULL;
  }
  // 1. Scaling
  const int scale = 256;
  unsigned char *scaledImage =
      reduceSize(colorImage, width, height, channels, scale);
  stbi_image_free(colorImage);

  // 2. Color space
  double *grayImage = grayscale(scaledImage, scale, scale, channels);
  free(scaledImage);

  // 3. Wavelet
  const int max_level = (int)(log(scale) * M_LOG2E + 1e-9);
  const int hash_level = (int)(log(sqrt(hashSize)) * M_LOG2E + 1e-9);
  const int level = max_level - hash_level;
  double *haar = createHaarImage(grayImage, scale, scale, level);
  free(grayImage);

  // 4. Truncation
  // const int t = 40;
  // double *truncated = truncateCoefficients(haar, scale, t);
  // free(truncated);

  // 5. Quantization
  // TODO: check hash size is power of 2
  int stride = (int)(sqrt(hashSize) + 1e-9);
  double *fingerprint = (double *)malloc(sizeof(double) * hashSize);
  for (int r = 0; r < stride; r++) {
    for (int c = 0; c < stride; c++) {
      fingerprint[(r) * (stride) + (c)] = haar[(r)*scale + (c)];
    }
  }
  unsigned char *hash = quantifyMedian(fingerprint, hashSize);
  free(fingerprint);

  // 6. Normalization

  return hash;
}

int main(int argc, char *argv[]) {
  if (argc < 3) {
    printf("Usage: %s <image1.jpg> <image2.jpg>\n", argv[0]);
    return 2;
  }

  int hashSize = 64;
  int err = 0;
  unsigned char *hash1 = getFingerprint(&err, argv[1], hashSize);
  if (err) {
    return err;
  }
  unsigned char *hash2 = getFingerprint(&err, argv[2], hashSize);
  if (err) {
    return err;
  }

  double diff = 0;
  for (int i = 0; i < hashSize; i++) {
    if (hash1[i] != hash2[i]) {
      diff += 1;
    }
  }
  printf("%.5f\n", (1 - (diff / hashSize)) * 100);

  free(hash1);
  free(hash2);

  return 0;
}
