#ifndef MD5_H
#define MD5_H

// Code adapted from http://sourceforge.net/projects/hashlib2plus

#include <string>
#include <stdio.h>
using namespace std;

// Define useful typedef, and padding array
typedef unsigned char *POINTER;

static unsigned char PADDING[64] = { 0x80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0 };

// Define constants for md5
#define S11 7
#define S12 12
#define S13 17
#define S14 22
#define S21 5
#define S22 9
#define S23 14
#define S24 20
#define S31 4
#define S32 11
#define S33 16
#define S34 23
#define S41 6
#define S42 10
#define S43 15
#define S44 21

// Define functions required for MD5
#define F(x, y, z) (((x) & (y)) | ((~x) & (z)))
#define G(x, y, z) (((x) & (z)) | ((y) & (~z)))
#define H(x, y, z) ((x) ^ (y) ^ (z))
#define I(x, y, z) ((y) ^ ((x) | (~z)))

#define ROTATE_LEFT(x, n) (((x) << (n)) | (( (unsigned int) x) >> (32-(n))))

#define FF(a, b, c, d, x, s, ac) { \
 (a) += F ((b), (c), (d)) + (x) + (unsigned long int)(ac); \
 (a) = ROTATE_LEFT ((a), (s)); \
 (a) += (b); \
  }
#define GG(a, b, c, d, x, s, ac) { \
 (a) += G ((b), (c), (d)) + (x) + (unsigned long int)(ac); \
 (a) = ROTATE_LEFT ((a), (s)); \
 (a) += (b); \
  }
#define HH(a, b, c, d, x, s, ac) { \
 (a) += H ((b), (c), (d)) + (x) + (unsigned long int)(ac); \
 (a) = ROTATE_LEFT ((a), (s)); \
 (a) += (b); \
  }
#define II(a, b, c, d, x, s, ac) { \
 (a) += I ((b), (c), (d)) + (x) + (unsigned long int)(ac); \
 (a) = ROTATE_LEFT ((a), (s)); \
 (a) += (b); \
  }

// Define struct to contain MD5 context
typedef struct {
	// state (ABCD)
	unsigned long int state[4];

	// number of bits, modulo 2^64 (lsb first)
	unsigned long int count[2];

	// input buffer
	unsigned char buffer[64];
} MD5_CTX;

class MD5 {
private:
	MD5_CTX ctx;

	string hashIt(void);
	string convToString(unsigned char *data);
	void updateContext(unsigned char *data, unsigned int len);
	void resetContext(void);

	void MD5Transform(unsigned long int state[4], unsigned char block[64]);
	void Encode(unsigned char* output, unsigned long int *input,
			unsigned int len);
	void Decode(unsigned long int *output, unsigned char *input,
			unsigned int len);
	void MD5_memcpy(POINTER output, POINTER input, unsigned int len);
	void MD5_memset(POINTER output, int value, unsigned int len);

public:
	MD5();
	~MD5();

	string getHashFromFile(string filename);

	void MD5Init(MD5_CTX* context);
	void MD5Update(MD5_CTX* context, unsigned char *input,
			unsigned int inputLen);
	void MD5Final(unsigned char digest[16], MD5_CTX* context);
};

#endif

