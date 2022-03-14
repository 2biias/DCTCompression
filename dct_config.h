#include <stdint.h>

/* Control the program using the four flags below */
#define DCT_M 4
#define DCT_L 8
#define FIXED_POINT_DATA 1
//#define PRINT_OUTPUT

#define DCT_N 512
#define COMPRESSIONRATIO (DCT_L / DCT_M)
#define DCT_ENCODED_N (DCT_N / COMPRESSIONRATIO)
#define DCT_WINDOWS (DCT_N / DCT_L)

#define FIXED_FRACTION 12
#define FIXED_INTEGER 4
#define FIXED_MULTIPLICATION(x, y) (((x) >> FIXED_FRACTION/2) * ((y) >> FIXED_FRACTION/2))

#define FLOAT2FIXED(x) ((short)((x) * (1 << FIXED_FRACTION)))
#define FIXED2FLOAT(x) (((float)(x)) / (1 << FIXED_FRACTION))

float determinant(int order, float dct_matrix_[order][order]);
void transpose(int order, float dct_matrix_[order][order]);
void cofactor(int order, float dct_matrix_[order][order]);
void calculate_DCT_matrix(float dct_matrix_[DCT_L][DCT_L]);
void calculate_DCT_matrix_reduced(uint8_t DCT_M_, uint8_t DCT_L_, float dct_matrix_[DCT_M_][DCT_L_]);
void calculate_DCT_inverse_matrix(int order, float dct_matrix_[order][order]);

#ifdef FIXED_POINT_DATA
static short encoded[DCT_ENCODED_N];
static short decoded[DCT_N];
void DCT_encode_data(short dct_matrix_[DCT_L][DCT_L], short samples_[static DCT_N], short encoded_[static DCT_ENCODED_N]);
void DCT_decode_data(short dct_inverse_matrix_[DCT_L][DCT_L], short decoded_[static DCT_N], short encoded_[static DCT_ENCODED_N]);
float RMSError(short samples_[DCT_N], short decoded_[DCT_N]);
#else
static float encoded[DCT_ENCODED_N];
static float decoded[DCT_N];
void DCT_encode_data(float dct_matrix_[DCT_L][DCT_L], float samples_[static DCT_N], float encoded_[static DCT_ENCODED_N]);
void DCT_decode_data(float dct_inverse_matrix_[DCT_L][DCT_L], float decoded_[static DCT_N], float encoded_[static DCT_ENCODED_N]);
float RMSError(float samples_[DCT_N], float decoded_[DCT_N]);
#endif

/* Misc */
void printFloat(float f);