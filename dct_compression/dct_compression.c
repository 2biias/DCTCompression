#include "contiki.h"
#include "sys/rtimer.h"
#include "dev/button-sensor.h"
#include "lib/sensors.h"
#include "math.h"
#include "mit200.h"
#include "dct_config.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define PI 3.141592653589793238462643383

float determinant(int order, float dct_matrix_[order][order]) {
    float a[order][order], ratio, det=1;
	int i, j, k;
    for(i=0; i < order; ++i) {
        for(j=0; j < order; ++j) {
            a[i][j] = dct_matrix_[i][j];
        }
    }
    /* Gauss elimination to get upper triangular matrix. */
    for(i=0; i < order; ++i)
    {
        if(a[i][i] == 0.0) {
            /*  
                Unfortunately Gauss elimination doesn't work with zeros in the diagonal
                even though the matrix is invertible. For now I assume the H matrix to
                not have zeros in the diagonal. Otherwise another technique could be
                implemented.
            */
            printf("Zero in diagonal - result is wrong\n");
            return 0;
        }
        for(j=i + 1; j < order; ++j) {
            ratio = a[j][i]/a[i][i];
            for(k=0; k < order; ++k) {
                a[j][k] = a[j][k] - ratio*a[i][k];
            }
        }
    }
    // Calculate the determinant
    for(i=0; i < order; ++i) {
         det = det * a[i][i];
    }
    return det;
}

void transpose(int order, float dct_matrix_[order][order]) {
    float temp;
    int i, j;
    for(i = 0; i < order; ++i) {
        for(j = 0; j < order; ++j) {
            if(i == j || i > j) continue;
            temp = dct_matrix_[i][j];
            dct_matrix_[i][j] = dct_matrix_[j][i];
            dct_matrix_[j][i] = temp;
        }
    }
}

void cofactor(int order, float dct_matrix_[order][order]) {
    int i, j, k, l;
    int suborder = order-1;
    float submatrix[suborder][suborder];
    
    /* Transpose the matrix */
    transpose(order, dct_matrix_);

    for(i = 0; i < order; ++i) {
        for(j = 0; j < order; ++j) {
            // Extract a minor
            int x = 0;
            int y = 0;
            for(k = 0; k < order; ++k) {
                if(k == i) continue;
                y = 0;
                for(l = 0; l < order; ++l) {
                    if(j == l) {
                        continue;
                    }
                    submatrix[x][y] = dct_matrix_[k][l];
                    y++;
                }
                x++;
            }
            // Find determinant of minor
            if( (i+j) % 2 )
                dct_matrix_[i][j] = (-1) * determinant(suborder, submatrix);
            else
                dct_matrix_[i][j] = determinant(suborder, submatrix);
        }
    }
}

void calculate_DCT_matrix(float dct_matrix_[DCT_L][DCT_L]) {
    uint8_t p, q;
    for(p = 0; p < DCT_L; ++p) {
        for(q = 0; q < DCT_L; ++q) {
            if(!p)
                dct_matrix_[p][q] = 1 / sqrt(DCT_L);
            else
                dct_matrix_[p][q] = sqrt(2/(float)(DCT_L))*cos((PI*(2*q+1)*p)/(2*DCT_L));
        }
    }
}

void calculate_DCT_matrix_reduced(uint8_t DCT_M_, uint8_t DCT_L_, float dct_matrix_[DCT_M_][DCT_L_]) {
    uint8_t p, q;
    for(p = 0; p < DCT_M_; ++p) {
        for(q = 0; q < DCT_L_; ++q) {
            if(!p)
                dct_matrix_[p][q] = 1 / sqrt(DCT_L_);
            else
                dct_matrix_[p][q] = sqrt(2/(float)(DCT_L))*cos((PI*(2*q+1)*p)/(2*DCT_L_));
        }
    }
}

void calculate_DCT_inverse_matrix(int order, float dct_matrix_[order][order]) {
    float det = determinant(order, dct_matrix_);
    if(det != 0) {
        cofactor(order, dct_matrix_);
        float factor = 1/det;
        int i, j;
        for(i = 0; i < order; ++i) {
            for(j = 0; j < order; ++j) {
                dct_matrix_[i][j] = factor * dct_matrix_[i][j];
            }
        }
    }
}

#ifdef FIXED_POINT_DATA

static short encoded[DCT_ENCODED_N] = {0};
static short decoded[DCT_N] = {0};

void DCT_encode_data(short dct_matrix_[DCT_L][DCT_L], short samples_[static DCT_N], short encoded_[static DCT_ENCODED_N]) {
    uint8_t i, j, k;
    for(i = 0; i < DCT_WINDOWS; ++i) {
        for(j = 0; j < DCT_M; ++j) {
            for(k = 0; k < DCT_L; ++k) {
                /* Use the hardware multiplier to do 16x16->32 multiplication (avoid type casting) */
                MACS = dct_matrix_[j][k];
                OP2 = samples_[i*DCT_L + k];
            }
            encoded_[i*DCT_M + j] = (RESHI << FIXED_INTEGER | RESLO >> FIXED_FRACTION);
            RESHI = 0;
            RESLO = 0;
        }
    }
}

void DCT_decode_data(short dct_inverse_matrix_[DCT_L][DCT_L], short decoded_[static DCT_N], short encoded_[static DCT_ENCODED_N]) {
    uint8_t i, j, k;
    for(i = 0; i < DCT_WINDOWS; ++i) {
        for(j = 0; j < DCT_L; ++j) {
            for(k = 0; k < DCT_M; ++k) {
                /* Use the hardware multiplier to do 16x16->32 multiplication (avoid type casting) */
                MACS = dct_inverse_matrix_[j][k];
                OP2 = encoded_[DCT_M*i + k];
            }
            decoded_[DCT_L*i + j] = (RESHI << FIXED_INTEGER | RESLO >> FIXED_FRACTION);
            RESHI = 0;
            RESLO = 0;
        }
    }
}

float RMSError(short samples_[DCT_N], short decoded_[DCT_N]) {
    uint16_t n;
    float sum = 0;
    float sample_float_, decoded_float_;

    for(n = 0; n < DCT_N; ++n) {
        sample_float_ = FIXED2FLOAT(samples_[n]);
        decoded_float_ = FIXED2FLOAT(decoded_[n]);
        sum += pow(sample_float_ - decoded_float_, 2);
    }
    return sqrt(sum / (float)(DCT_N));
}

#else // FIXED_POINT_DATA

static float encoded[DCT_ENCODED_N] = {0};
static float decoded[DCT_N] = {0};

void DCT_encode_data(float dct_matrix_[DCT_L][DCT_L], float samples_[static DCT_N], float encoded_[static DCT_ENCODED_N]) {
    uint16_t i, j, k;
    for(i = 0; i < DCT_WINDOWS; ++i) {
        for(j = 0; j < DCT_M; ++j) {
            for(k = 0; k < DCT_L; ++k) {
                encoded_[i*DCT_M + j] += dct_matrix_[j][k] * samples_[i*DCT_L + k];
            }
        }
    }
}

void DCT_decode_data(float dct_inverse_matrix_[DCT_L][DCT_L], float decoded_[static DCT_N], float encoded_[static DCT_ENCODED_N]) {
    uint16_t i, j, k;
    for(i = 0; i < DCT_WINDOWS; ++i) {
        for(j = 0; j < DCT_L; ++j) {
            for(k = 0; k < DCT_M; ++k) {
                decoded_[DCT_L*i + j] += dct_inverse_matrix_[j][k] * encoded_[DCT_M*i + k];
            }
        }
    }
}

float RMSError(float samples_[DCT_N], float decoded_[DCT_N]) {
    uint16_t n;
    float sum = 0;
    for(n = 0; n < DCT_N; ++n) {
        sum += pow(samples_[n] - decoded_[n], 2);
    }
    return sqrt(sum / (float)(DCT_N));
}

#endif // FIXED_POINT_DATA


/* Prints float WITH sign */
void printFloat(float f) {
    char *format;
    int v = f * 10000;
    if (v < 0) {
        format = "-%d.%04d";
        v = -v;
    } else {
        format = "%d.%04d";
    }
    printf(format, v/10000, v%10000);
}

/*-----------------------------------------------------------*/
PROCESS(encode_process, "DCT encoder process");

AUTOSTART_PROCESSES(&encode_process);
/*-----------------------------------------------------------*/

// Define global variables
PROCESS_THREAD(encode_process, ev, data)
{
    int y, x;
    rtimer_clock_t t1, t2;
    unsigned long time_ms;
    SENSORS_ACTIVATE(button_sensor);
    PROCESS_BEGIN();
    while(1) {
        PROCESS_WAIT_EVENT_UNTIL(ev == sensors_event && data == &button_sensor);

        // Watchdog is disabled/enabled because the floating point version takes > 1000ms
        watchdog_stop();

        float (*dct_matrix)[DCT_L] = malloc(sizeof(float[DCT_L][DCT_L]));
        float (*dct_inverse_matrix)[DCT_L] = malloc(sizeof(float[DCT_L][DCT_L]));
        
        /* Calulate H matrix */
        t1 = RTIMER_NOW();
        calculate_DCT_matrix(dct_matrix);
        t2 = RTIMER_NOW();
        t2 = t2 - t1;
        time_ms = ((unsigned long)t2 * 1000) / RTIMER_SECOND;

        printf("Time spend calc. H-matrix: %lums\n", time_ms);
        for(x = 0; x < DCT_L; ++x) {
            for(y = 0; y < DCT_L; ++y) {
                dct_inverse_matrix[x][y] = dct_matrix[x][y];
            }
        }
        calculate_DCT_inverse_matrix(DCT_L, dct_inverse_matrix);

        #ifdef FIXED_POINT_DATA
        short (*dct_matrix_fixed)[DCT_L] = malloc(sizeof(short[DCT_L][DCT_L]));
        for(x = 0; x < DCT_L; ++x) {
            for(y = 0; y < DCT_L; ++y) {
                dct_matrix_fixed[x][y] = FLOAT2FIXED(dct_matrix[x][y]);
            }
        }
        free(dct_matrix);
        short (*dct_inverse_matrix_fixed)[DCT_L] = malloc(sizeof(short[DCT_L][DCT_L]));
        for(x = 0; x < DCT_L; ++x) {
            for(y = 0; y < DCT_L; ++y) {
                dct_inverse_matrix_fixed[x][y] = FLOAT2FIXED(dct_inverse_matrix[x][y]);
            }
        }
        free(dct_inverse_matrix);
        #endif

        /* Encode data */
        t1 = RTIMER_NOW();
        #ifdef FIXED_POINT_DATA
        DCT_encode_data(dct_matrix_fixed, mit200, encoded);
        #else
        DCT_encode_data(dct_matrix, mit200, encoded);
        #endif
        t2 = RTIMER_NOW();
        t2 = t2 - t1;
        time_ms = ((unsigned long)t2 * 1000) / RTIMER_SECOND;
        printf("Time spend encoding: %lums\n", time_ms);

        /* Decode data */
        t1 = RTIMER_NOW();
        #ifdef FIXED_POINT_DATA
        DCT_decode_data(dct_inverse_matrix_fixed, decoded, encoded);
        #else
        DCT_decode_data(dct_inverse_matrix, decoded, encoded);
        #endif
        t2 = RTIMER_NOW();
        t2 = t2 - t1;
        time_ms = ((unsigned long)t2 * 1000) / RTIMER_SECOND;
        printf("Time spend decoding: %lums\n", time_ms);

        printf("Raw data:\n");
        for(y = 0; y < 32; ++y) {
            #ifdef FIXED_POINT_DATA
            printFloat(FIXED2FLOAT(mit200[y]));
            #else
            printFloat(mit200[y]);
            #endif
            putchar('\n');
        }

        printf("Encoded data:\n");
        for(y = 0; y < 16; ++y) {
            #ifdef FIXED_POINT_DATA
            printFloat(FIXED2FLOAT(encoded[y]));
            #else
            printFloat(encoded[y]);
            #endif
            putchar('\n');
        }

        printf("Decoded data:\n");
        for(y = 0; y < 32; ++y) {
            #ifdef FIXED_POINT_DATA
            printFloat(FIXED2FLOAT(decoded[y]));
            #else
            printFloat(decoded[y]);
            #endif
            putchar('\n');
        }

        printf("RMS Error: ");
        printFloat(RMSError(mit200, decoded));
        putchar('\n');

        #ifdef PRINT_OUTPUT
        putchar('[');
        for(y = 0; y < DCT_N; ++y) {
            #ifdef FIXED_POINT_DATA
                printFloat(FIXED2FLOAT(decoded[y]));
            #else
                printFloat(decoded[y]);
            #endif
            if(y != DCT_N-1) putchar(' ');
        }
        printf("]\n");
        #endif

        for(y = 0; y < DCT_N; ++y) decoded[y] = 0;
        for(y = 0; y < DCT_ENCODED_N; ++y) encoded[y] = 0;

        #ifdef FIXED_POINT_DATA
        free(dct_matrix_fixed);
        free(dct_inverse_matrix_fixed);
        #else
        free(dct_matrix);
        free(dct_inverse_matrix);
        #endif

        watchdog_start();
    }
    PROCESS_END();
}
